#!/bin/bash -f
# -f flag unsets wildcard expansion which is needed for the move operations
USAGE=" 

Usage: $0 [-nfh] [-c halorc] projectrc

Positional arguments
  
  projectrc  - project configuration file

Options

  -n         - dry run
  -f         - force execution of all steps
  -c         - halo runtime configuration. Default: ~/.halorc

"

# Getting the path to the current script
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESYNCMATES=$DIR/resyncMates.pl

# Force move operations
MV="mv -f"

# Setup environment
DRY_RUN=false
FORCE=false

E_PARAM_ERR=250       # If less than 2 params passed to function.
P_RUN=true            # Run command
P_NORUN=false         # Skip command
STAGING_DIR=tx        # Staging directory in which to generate temporary files
JAVA_OPTS=-Xmx3g
N_CORES=1
# Initialize stdout and stderr variables; can be overloaded in config
# file. These variables set where *program output* logging goes;
# pipeline status commands will still be printed to console
LOGFILE=/dev/stdout
ERRFILE=/dev/stderr

# Halo configuration file
halorc=~/.halorc
# Check user input
while getopts :nfhc: flag; do
  case $flag in
    n)
      dry_run=true;
      ;;
    f)
      force=true;
      ;;
    c) 
      halorc=$OPTARG;
      ;;
    h)
      echo "$USAGE";
      exit;
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "$USAGE";
      exit;
      ;;
  esac
done
shift $(( OPTIND - 1 ));

# Exit if no input
if [ $# -ne 1 ]; then
    echo "ERROR: Wrong number of arguments"
    echo "$USAGE"
    exit;
fi

# Functions
# run_command: emulate dry run behaviour
# NOTE: currently not used!
run_command () {
    if [ -z "$1" ]
    then
	echo $(date) "No command passed to run_command"
	exit
    fi
    command=$1
    output=$2
    if [ "$dry_run" ]; then
	echo $(date) "Command: " $command >> $LOGFILE
    else
	up_to_date $output
	retval=$?
	if [ ! "$retval" ]; then
	    echo $(date) "command up to date; skipping"
	    echo -e "\t" $command >> $LOGFILE
	else
	    $command
	fi
    fi
}

# up_to_date: check if job is up to date. Very simple rules: if only
# output file given, we assume that job is done if file exists. If
# pre_output file given, compare time stamps and return 0 if output
# older than pre_output (emulates make).
up_to_date () {
    if [ $force ]; then
	echo "Force flag set; running analysis" >> $LOGFILE
	return 0
    fi
    if [ -z "$1" ]
    then
	echo "No parameter passed to up_to_date" >> $ERRFILE
	exit
    fi
    up_to_date_output=$1
    if [ ! -e $up_to_date_output ]; then
	echo "up_to_date: No such file $up_to_date_output; running analysis" >> $ERRFILE
	return 0
    fi
    if [ ! -z "$2" ]; then
	up_to_date_pre_output=$2
	if [ -e "$up_to_date_pre_output" ]; then
	    up_to_date_output_time=`stat -c %Y $up_to_date_output`
	    up_to_date_pre_output_time=`stat -c %Y $up_to_date_pre_output`
	    if [ $up_to_date_output_time -gt $up_to_date_pre_output_time ];
	    then
		return 1
	    else
		return 0
	    fi
	fi
    fi
    return 1
}

##################################################
# Start processing samples
##################################################
# Source the project runtime configuration. 
projectrc=$1
if [ ! -e $projectrc ]; then
    echo $(date) "no project configuration file found; no such file $projectrc" >> $ERRFILE
    exit
fi
source $projectrc

# Source the halo runtime configuration.
if [ ! -e $halorc ]; then
    echo $(date) "no halo configuration file found; no such file $halorc" >> $ERRFILE
    exit
fi
source $halorc
# Set parallel options
PARALLEL_OPTS="-j $N_CORES"

if [ ! "$samples" ]; then
    echo  $(date) "No samples provided; please provide either sample names or a file name listing sample names in the project configuration file" >> $ERRFILE
    echo "$USAGE" >> $ERRFILE
    exit
fi

if [ ! "$indir" ]; then
    echo  $(date) "No indir provided; please provide indir in the project configuration file" >> $ERRFILE
    echo "$USAGE" >> $ERRFILE
    exit
fi
# Init logging
echo $(date) > $LOGFILE
echo $(date) > $ERRFILE

# Make sure indir is absolute path
if [ "`echo $indir | cut -c1`" != "/" ]; then
    echo  $(date) "$indir is not an absolute path; please use absolute path names for input directory" >> $ERRFILE
    exit
fi

# Set input directories
if [ -f "$samples" ]; then
    echo  $(date) "Reading input directories from file $samples"
    samples=`cat $samples`
else
    echo $(date) "Assuming input directories passed as argument list"
fi
sample_regexp=`echo $samples | sed -e s"/ /\|/"g`

# Find input files based on regular expressions that include sample names
# Here assuming casava-based file names

echo "$(date) Finding files with command 'find $indir -regextype posix-extended -regex \".*(${sample_regexp}).*${READ1_REGEXP}.fastq.gz\"'"
infiles=`find $indir -regextype posix-extended -regex ".*(${sample_regexp}).*${READ1_REGEXP}.fastq.gz"`

echo $(date) "Going to run pipeline on "
for f in $infiles; do 
    echo -e "\t" ${f%.fastq.gz} 
done
read1=$infiles
read2=`echo $infiles | sed -e "s/${READ1_REGEXP}/${READ2_REGEXP}/g"`

# Setup staging directories
for f in $read1 $read2; do
    TMPDIR=`dirname $f`/tx
    mkdir -p $TMPDIR
done

##################################################
# Pipeline code
# The various steps follow a presentation on HaloPlex
# analysis obtained from Agilent
##################################################

##############################
# Primary analysis
##############################
# 1. QC
command=""
for f in $read1 $read2; do
    outdir=`dirname $f`/fastqc
    outfile=$outdir/`basename ${f%.fastq*}`_fastqc/summary.txt
    up_to_date  $outfile $f
    if [ $? = 1 ]; then continue; fi
    echo "running fastqc"  >> $LOGFILE
    mkdir -p $outdir
    cmd="$FASTQC $f -o ${outdir}"
    command="$command\n$cmd"
done
echo -e $(date) 1. QC section 
echo -e $(date) $command
echo -e "$command" | $PARALLEL $PARALLEL_OPTS  >> $LOGFILE 2>> $ERRFILE

# 2a. Trim adapter sequence
command=""
trimfiles=""
for f in $read1; do
    trimfiles="$trimfiles ${f%.fastq.gz}.trimmed.fastq.gz"
    TMP=`dirname $f`/$STAGING_DIR/`basename $f`
    up_to_date ${f%.fastq.gz}.trimmed.fastq.gz $f
    if [ $? = 1 ]; then continue; fi
    cmd="$CUTADAPT $CUTADAPT_OPTS -a $THREEPRIME $f -o ${TMP%.fastq.gz}.trimmed.fastq.gz > ${f%.fastq.gz}.trimmed.fastq.cutadapt_metrics && $MV ${TMP%.fastq.gz}.trimmed.fastq.gz ${f%.fastq.gz}.trimmed.fastq.gz"
    command="$command\n$cmd"
done
for f in $read2; do
    trimfiles="$trimfiles ${f%.fastq.gz}.trimmed.fastq.gz"
    TMP=`dirname $f`/$STAGING_DIR/`basename $f`
    up_to_date ${f%.fastq.gz}.trimmed.fastq.gz $f
    if [ $? = 1 ]; then continue; fi
    cmd="$CUTADAPT $CUTADAPT_OPTS -a $FIVEPRIME $f -o ${TMP%.fastq.gz}.trimmed.fastq.gz > ${f%.fastq.gz}.trimmed.fastq.cutadapt_metrics && $MV ${TMP%.fastq.gz}.trimmed.fastq.gz ${f%.fastq.gz}.trimmed.fastq.gz"
    command="$command\n$cmd"
done
echo -e $(date) 2a. Adapter trimming
echo -e $(date) $command
echo -e "$command" | $PARALLEL $PARALLEL_OPTS  >> $LOGFILE 2>> $ERRFILE

# 2b. Resync mates - sometimes cutadapt cuts reads down to 0, so there
# are some reads without mates
# NB: this process is memory hungry so currently set max_proc = 1
sample_pfx=`for f in $read1; do echo ${f%_${READ1_REGEXP}.fastq.gz}; done`
RESYNC_PROC=$N_CORES
RESYNC_MAX_PROC=1
if [ "$RESYNC_PROC" -gt "$RESYNC_MAX_PROC" ]; then RESYNC_PROC=$RESYNC_MAX_PROC; fi
command=""
syncfiles=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    TMP=$OUTDIR/$STAGING_DIR/`basename $f`
    syncfiles="$syncfiles ${TMP}_${READ1_REGEXP}.trimmed.sync.fastq.gz ${TMP}_${READ2_REGEXP}.trimmed.sync.fastq.gz"
    up_to_date ${f}_${READ1_REGEXP}.trimmed.sync.fastq.gz ${f}_${READ1_REGEXP}.trimmed.fastq.gz
    if [ $? = 1 ]; then continue; fi
    echo $(date) resyncing reads for ${f}_${READ1_REGEXP}.trimmed.fastq.gz, ${f}_${READ2_REGEXP}.trimmed.fastq.gz
    cmd="$RESYNCMATES -i ${f}_${READ1_REGEXP}.trimmed.fastq.gz -j ${f}_${READ2_REGEXP}.trimmed.fastq.gz -o ${TMP}_${READ1_REGEXP}.trimmed.sync.fastq.gz -p ${TMP}_${READ2_REGEXP}.trimmed.sync.fastq.gz && $MV ${TMP}_${READ1_REGEXP}.trimmed.sync.fastq.gz $OUTDIR && $MV ${TMP}_${READ2_REGEXP}.trimmed.sync.fastq.gz $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) 2b. Resync mates
echo -e $(date) $command
echo -e "$command" | $PARALLEL -j $RESYNC_PROC  >> $LOGFILE 2>> $ERRFILE

##############################
# Mapping - secondary analysis
##############################

# 3. Align sequences with bwa. Here we run command sequentially since
# bwa takes care of parallelization. From now on we run at sample level.
# Remove staging directory from syncfiles
MAX_PROC=`cat /proc/cpuinfo | grep processor | tail -1 | awk '{print $3}'`
syncfiles=`echo $syncfiles | sed -e "s/$STAGING_DIR\///g;"`
echo -e $(date) 3. Alignment
for f in $syncfiles; do
    OUTDIR=`dirname $f`
    TMP=$OUTDIR/$STAGING_DIR/`basename $f`
    up_to_date ${f%.fastq.gz}.sai $f
    if [ $? = 1 ]; then continue; fi
    echo $(date) aligning reads $f
    echo $(date) "$BWA aln -t $MAX_PROC $BWA_REF $f > ${TMP%.fastq.gz}.sai" >> $ERRFILE
    $BWA aln -t $MAX_PROC $BWA_REF $f > ${TMP%.fastq.gz}.sai 2>> $ERRFILE && $MV ${TMP%.fastq.gz}.sai $OUTDIR
done

# 4. Pair reads
# http://bio-bwa.sourceforge.net/bwa.shtml: "...the sampe command uses ~5.4GB"
SAMPE_PROC=$N_CORES
SAMPE_MAX_PROC=`cat /proc/meminfo | grep MemTotal | awk  '{print int($2/1e6/5.4)}'`
if [ "$N_CORES" -gt "$SAMPE_MAX_PROC" ]; then SAMPE_PROC=$SAMPE_MAX_PROC; fi
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    TMP=$OUTDIR/$STAGING_DIR/`basename $f`
    label=`basename $f`
    extension="trimmed.sync"
    up_to_date $f.sam ${f}_${READ1_REGEXP}.${extension}.sai
    if [ $? = 1 ]; then continue; fi
    echo $(date) pairing reads for sample $f
    cmd="$BWA sampe -A -P -r \"@RG\tID:${label}\tSM:${label}\tPL:Illumina\tCN:Agilent\" $BWA_REF ${f}_${READ1_REGEXP}.${extension}.sai ${f}_${READ2_REGEXP}.${extension}.sai ${f}_${READ1_REGEXP}.${extension}.fastq.gz ${f}_${READ2_REGEXP}.${extension}.fastq.gz > ${TMP}.sam && $MV ${TMP}.sam $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) 4. Pair reads
echo -e $(date) $command
echo -e $command | $PARALLEL -j $SAMPE_PROC >> $LOGFILE 2>> $ERRFILE

# 5. Generate bam file
# Records in RAM = 250000 times # Gb RAM - cf http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:___A_Picard_program_that_sorts_its_output_SAM.2FBAM_file_is_taking_a_very_long_time_and.2For_running_out_of_memory.__What_can_I_do.3F
# Similar to SAMPE we seem to run out of resources - in this case CPU
SORT_PROC=$SAMPE_PROC
SORTSAM_OPTS="SO=coordinate MAX_RECORDS_IN_RAM=750000"
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    TMP=$OUTDIR/$STAGING_DIR/`basename $f`
    JAVA_EXTRA_OPTS=-Djava.io.tmpdir=$OUTDIR/tmp
    up_to_date $f.sort.bam $f.bam
    if [ $? = 1 ]; then continue; fi
    echo "$SAMTOOLS view -bS $f.sam > ${TMP}.bam && java $JAVA_EXTRA_OPTS $JAVA_OPTS -jar $PICARD_HOME/SortSam.jar $SORTSAM_OPTS INPUT=${TMP}.bam OUTPUT=${TMP}.sort.bam && echo \"$f.bam removed to save space; see $f.sort.bam\" > ${TMP}.bam && $MV ${TMP}*ba* $OUTDIR && echo \"$f.sam removed to save space; see $f.sort.bam\" > $f.sam && samtools index $f.sort.bam"
    cmd="$SAMTOOLS view -bS $f.sam > ${TMP}.bam && java  $JAVA_EXTRA_OPTS $JAVA_OPTS -jar $PICARD_HOME/SortSam.jar $SORTSAM_OPTS INPUT=${TMP}.bam OUTPUT=${TMP}.sort.bam && echo \"$f.bam removed to save space; see $f.sort.bam\" > ${TMP}.bam && $MV ${TMP}*ba* $OUTDIR && echo \"$f.sam removed to save space; see $f.sort.bam\" > $f.sam && samtools index $f.sort.bam"
    command="$command\n$cmd"
done
echo -e $(date) 5. Generate sorted bam file
echo -e $(date) $command
echo -e $command | $PARALLEL -j $SORT_PROC >> $LOGFILE 2>> $ERRFILE

# 6. Generate various metrics for the bamfiles
# No staging directory is used here
command=""
for f in $sample_pfx; do
    input=$f.sort.bam
    # Alignment metrics
    up_to_date ${input%.bam}.align_metrics $input
    if [ $? = 0 ]; then 
	echo $(date) generating alignment metrics for $input
	cmd="java -jar $PICARD_HOME/CollectAlignmentSummaryMetrics.jar INPUT=$input OUTPUT=${input%.bam}.align_metrics REFERENCE_SEQUENCE=$REF"
	command="$command\n$cmd"
    fi
    # Insert size metrics
    up_to_date ${input%.bam}.insert_metrics $input
    if [ $? = 0 ]; then 
	echo $(date) generating insert size metrics for $input
	cmd="java -jar $PICARD_HOME/CollectInsertSizeMetrics.jar INPUT=$input OUTPUT=${input%.bam}.insert_metrics HISTOGRAM_FILE=${input%.bam}.insert_hist REFERENCE_SEQUENCE=$REF"
	command="$command\n$cmd"
    fi
    # Duplication statistics
    up_to_date ${input%.bam}.dup_metrics $input
    if [ $? = 0 ]; then
	echo $(date) generating duplication metrics for $input
	cmd="java -jar $PICARD_HOME/MarkDuplicates.jar INPUT=$input METRICS_FILE=${input%.bam}.dup_metrics OUTPUT=${input%.bam}.dup.bam"
	command="$command\n$cmd"
    fi
    # Hybrid selection metrics
    up_to_date ${input%.bam}.hs_metrics
    if [ $? = 0 ]; then 
	if [ ! $BAIT_INTERVALS_FILE ] || [ ! $TARGET_INTERVALS_FILE ]; then
	    echo $(date) "Bait/target file missing; skipping hybrid metrics calculation for $f"
	else 
	    echo $(date) generating hybrid selection metrics for $input
	    cmd="java -jar $PICARD_HOME/CalculateHsMetrics.jar INPUT=$input OUTPUT=${input%.bam}.hs_metrics BAIT_INTERVALS=$BAIT_INTERVALS_FILE TARGET_INTERVALS=$TARGET_INTERVALS_FILE REFERENCE_SEQUENCE=$REF"
	    command="$command\n$cmd"
	fi
    fi
    # Samstat statistics
    up_to_date $input.html $input
    if [ $? = 0 ]; then
	echo $(date) generating samstat metrics for $input
	cmd="$SAMSTAT $input"
	command="$command\n$cmd"
    fi
done
echo -e $(date) 6. Calculate metrics
echo -e $(date) $command
echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE
##############################
# Variant calling, tertiary analysis
##############################

# 7. Raw variant calling 
# -L: interval_list(?) type input - bed format doesn't seem to work
UNIFIED_GENOTYPER_OPTS="-T UnifiedGenotyper --dbsnp $DBSNP -stand_call_conf 30.0 -stand_emit_conf 10.0  --downsample_to_coverage 30 --output_mode EMIT_VARIANTS_ONLY -glm BOTH -nt $N_CORES -R $REF -L $TARGET_INTERVALS_FILE"
echo -e $(date) 7. Raw variant calling
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.bam
    TMP=$OUTDIR/$STAGING_DIR/`basename $input`
    output=${TMP%.bam}.BOTH.raw.vcf
    up_to_date ${input%.bam}.BOTH.raw.vcf $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) generating raw variant calls for $f
    samtools index $f.sort.bam && java -jar $GATK $UNIFIED_GENOTYPER_OPTS -I $input -o $output && $MV $output $OUTDIR >> $LOGFILE 2>> $ERRFILE
done

# 8. Generate realignment intervals 
DEEP_COVERAGE_OPTS="--mismatchFraction 0.30 --maxIntervalSize 650"
REALIGNMENT_TARGET_CREATOR_OPTS="-T RealignerTargetCreator -L $TARGET_INTERVALS_FILE -known $MILLS -known $THOUSANDG_OMNI -nt $N_CORES -R $REF $DEEP_COVERAGE_OPTS"
echo -e $(date) 8. Realignment - generate realign target intervals
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.intervals $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) generating realignment intervals for $input
    java -jar $GATK $REALIGNMENT_TARGET_CREATOR_OPTS -I $input -known ${input%.bam}.BOTH.raw.vcf -o ${output%.bam}.intervals && $MV ${output%.bam}.intervals $OUTDIR >> $LOGFILE 2>> $ERRFILE
done

# 9. Indel realignment 
#
# According to the documentation, IndelRealigner does *not* support
# option -nt (number of threads). Nevertheless, this process cranks up
# the CPU to as many cores as are available. I am missing some option
# here.

# NB: Could be due to Javas GC: http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_Why_does_a_Picard_program_use_so_many_threads.3F

# One or several of DEEP_COVERAGE_OPTS seems to screw up the output
DEEP_COVERAGE_OPTS="--maxReadsInMemory 300000 --maxReadsForRealignment 500000 --maxReadsForConsensuses 500 --maxConsensuses 100"
INDEL_REALIGNER_OPTS="-T IndelRealigner -known $MILLS -known $THOUSANDG_OMNI  -R $REF"
command=""
echo -e $(date) 9. Indel realignment
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.realign.bam $input
    if [ $? = 1 ]; then continue; fi
    echo -e $(date) Realigning round indels for $f
    java -jar $GATK $INDEL_REALIGNER_OPTS -I $input --targetIntervals ${input%.bam}.intervals -known ${input%.bam}.BOTH.raw.vcf -o ${output%.bam}.realign.bam && $MV "${output%.bam}.realign.ba*" $OUTDIR >> $LOGFILE 2>> $ERRFILE
    # command="$command\n$cmd"
done

# echo -e $(date) $command
# echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE

# 10. Base recalibration
# NB: currently BaseRecalibrator does *not* support multiple threads
RECALIBRATOR_OPTS="-T BaseRecalibrator -R $REF -L $TARGET_INTERVALS_FILE --knownSites $DBSNP --knownSites $MILLS --knownSites $THOUSANDG_OMNI"
echo -e $(date) 10. Base recalibration
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.recal_data.grp $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) realigning $input
    cmd="java -jar $GATK $RECALIBRATOR_OPTS -I $input -o ${output%.bam}.recal_data.grp && $MV ${output%.bam}.recal_data.grp $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE

# 11. Recalculate base quality score
PRINT_READS_OPTS="-T PrintReads -R $REF"
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.recal.bam $input
    if [ $? = 1 ]; then continue; fi
    cmd="java -jar $GATK $PRINT_READS_OPTS -I $input -BQSR ${input%.bam}.recal_data.grp -o ${output%.bam}.recal.bam && $MV "${output%.bam}.recal.ba*" $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) 11. Recalculate base quality score
echo -e $(date) $command
echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE


# 12. Clip 5 bp in 5' on all reads. These bases are reference bias due
# to the restriction enzyme site. 
CLIP_READS_OPTS="-T ClipReads --cyclesToTrim 1-5 --clipRepresentation WRITE_NS -R $REF"
command=""
echo -e $(date) 12. Clip 5 bp from 5 prime end of reads
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.recal.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.clip.bam $input
    if [ $? = 1 ]; then continue; fi
    cmd="java -jar $GATK $CLIP_READS_OPTS -I $input -o ${output%.bam}.clip.bam && $MV ${output%.bam}.clip.ba* $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE

# 13. Make final calls with GATK UnifiedGenotyper
#     FIXME: change to HaplotypeCaller when appropriate
UNIFIEDGENOTYPER_OPTS="-T UnifiedGenotyper --dbsnp $DBSNP -filterMBQ -A AlleleBalance -stand_call_conf 30.0 -stand_emit_conf 10.0 -dt NONE --output_mode EMIT_VARIANTS_ONLY -glm BOTH -nt $N_CORES -L $TARGET_INTERVALS_FILE"
echo -e $(date) 13. Final variant calling with UnifiedGenotyper
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.recal.clip.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.BOTH.final.vcf $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) generating final variant calls for $f with UnifiedGenotyper
    java -jar $GATK $UNIFIED_GENOTYPER_OPTS -I $input -o ${output%.bam}.BOTH.final.vcf && $MV "${output%.bam}.BOTH.final.vcf*" $OUTDIR >> $LOGFILE 2>> $ERRFILE
done

# 14. Making final calls with samtools
SAMTOOLS_OPTS="-A -u -f $REF"
BCFTOOLS_OPTS="-v -c"
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.recal.clip.bam
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.bam}.BOTH.final.samtools.vcf $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) generating final variant calls for $f with samtools
    cmd="samtools mpileup $SAMTOOLS_OPTS $input > ${output%.bam}.BOTH.final.samtools.mpileup && bcftools view $BCFTOOLS_OPTS -g ${output%.bam}.BOTH.final.samtools.mpileup > ${output%.bam}.BOTH.final.samtools.vcf && $MV ${output%.bam}.BOTH.final.samtools.vcf $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) 14. Final variant calling with samtools 
echo -e $(date) $command
echo -e $command | $PARALLEL $PARALLEL_OPTS >> $LOGFILE 2>> $ERRFILE

# 15. Run final variant filtration
#
# MQ Root mean square of the mapping quality of reads covering (measure of varying quantity).
# MQ0 Number of reads overlapping with zero quality.
# DP Total (unfiltered) depth for called position.
# QUAL The Phred scale probability of being REF if ALT is called (and ALT if REF is called).
# QD QUAL / DP (variant confidence divided by unfiltered depth. Low scores are indications of false positives).

VARIANTFILTRATION_OPTS="-T VariantFiltration --clusterWindowSize 10 --clusterSize 3 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < 10\" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" --filterExpression \"QUAL > 30.0 && QUAL < 50.0\" --filterName \"LowQual\" --filterExpression \"QD < 1.5\" --filterName \"LowQD\" -R $REF "

command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.recal.clip.BOTH.final.vcf
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.vcf}.filtered.vcf $input
    if [ $? = 1 ]; then continue; fi
    cmd="java -jar $GATK $VARIANTFILTRATION_OPTS --variant $input -o ${output%.vcf}.filtered.vcf && $MV "${output%.vcf}.filtered.vcf*" $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) 15. Perform last variant filtration
echo -e $(date) "$command"
echo -e "$command" | $PARALLEL $PARALLEL_OPTS  >> $LOGFILE 2>> $ERRFILE

# 16. Run variant evaluation
# FIXME: should add  target_intervals
VARIANTEVAL_OPTS="-T VariantEval -R $REF --dbsnp $DBSNP -ST Filter -l INFO --doNotUseAllStandardModules --evalModule CompOverlap --evalModule CountVariants --evalModule GenotypeConcordance --evalModule TiTvVariantEvaluator --evalModule ValidationReport --stratificationModule Filter "
echo -e $(date) 16. Run variant evaluation
command=""
for f in $sample_pfx; do
    OUTDIR=`dirname $f`
    input=$f.sort.realign.recal.clip.BOTH.final.filtered.vcf
    output=$OUTDIR/$STAGING_DIR/`basename $input`
    up_to_date ${input%.vcf}.eval_metrics $input
    if [ $? = 1 ]; then continue; fi
    echo $(date) "Evaluate variants for $f"
    cmd="java -jar $GATK $VARIANTEVAL_OPTS --eval $input -o ${output%.vcf}.eval_metrics && $MV ${output%.vcf}.eval_metrics $OUTDIR"
    command="$command\n$cmd"
done
echo -e $(date) "$command"
echo -e "$command" | $PARALLEL $PARALLEL_OPTS  >> $LOGFILE 2>> $ERRFILE
