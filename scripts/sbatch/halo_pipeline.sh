#!/bin/bash

USAGE="Usage: $0 [-n] indir samples"

# Parameters that govern general behaviour
READ1_REGEXP="R1_001"    # Regexp used for 1st reads
READ2_REGEXP="R2_001"    # Regexp used for paired sequences
N_CORES=5
LOGFILE="halo_pipeline.out"
ERRFILE="halo_pipeline.err"

echo $(date) > $LOGFILE
echo $(date) > $ERRFILE

# Software config
# Modify for specific versions
PARALLEL=/usr/bin/parallel
PIGZ=/usr/bin/pigz
FASTQC=fastqc
CUTADAPT=cutadapt
CUTADAPT_OPTS="-m 50"

# Alignment options
BWA=bwa
BWA_HG19=/datad/biodata/genomes/Hsapiens/hg19/bwa/hg19.fa
BWA_REF=$BWA_HG19
REF=/datad/biodata/genomes/Hsapiens/hg19/hg19.fa

# Samtools
SAMTOOLS=samtools

# Adapter sequences. These correspond to TruSeq adapter sequence.
# THREEPRIME is found in the three-prime end of read 1, FIVEPRIME
# revcomp in the end of read 2
THREEPRIME="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
FIVEPRIME="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

# Setup environment
DRY_RUN=false
FORCE=false

E_PARAM_ERR=250    # If less than 2 params passed to function.
P_RUN=true            # Run command
P_NORUN=false         # Skip command

# Check user input
while getopts n:f flag; do
  case $flag in
    n)
      dry_run=true;
      ;;
    f)
      force=true;
      ;;
    ?)
      echo $usage
      exit;
      ;;
  esac
done
shift $(( OPTIND - 1 ));

# Exit if no input
if [ $# -eq 0 ]; then
    echo $usage
    exit;
fi

# Functions
# run_command: emulate dry run behaviour
run_command () {
    if [ -z "$1" ]
    then
	echo $(date) "No command passed to run_command"
	exit
    fi
    command=$1
    output=$2
    echo $output
    if [ "$dry_run" ]; then
	echo $(date) "Command: " $command
    else
	up_to_date $output
	retval=$?
	if [ ! "$retval" ]; then
	    echo $(date) "command up to date; skipping"
	    echo -e "\t" $command
	else
	    $command
	fi
    fi
}

# up_to_date: check if job is up to date. Very simple rule: if output file exists we assume that job is done, unless force flag is passed
up_to_date () {
    if [ -z "$1" ]
    then
	echo "No parameter passed to up_to_date"
	exit
    fi
    output=$1
    echo $(date) checking for file $output
    if [ -e "$output" ];
    then
	return 1
    else
	return 0
    fi
}

##################################################
# Start processing samples
##################################################
indir=$1
samples="${@:2}"
if [ ! -d "$indir" ]; then
    echo  $(date) "$indir not a directory"
    echo $USAGE
    exit
fi
if [ ! "$samples" ]; then
    echo  $(date) "No samples provided; please provide either sample names or a file listing sample names"
    echo $USAGE
    exit
fi
# Make sure indir is absolute path
if [ "`echo $indir | cut -c1`" != "/" ]; then
    echo  $(date) "$indir is not an absolute path; please use absolute path names for input directory"
    exit
fi

# Set input directories
if [ -f "$samples" ]; then
    echo  $(date) "Reading input directories from file $args"
    samples=`cat $samples`
else
    echo $(date) "Assuming input directories passed as argument list"
fi
sample_regexp=`echo $samples | sed -e s"/ /\|/"g`

# Find input files based on regular expressions that include sample names
# Here assuming casava-based file names
echo $(date) Finding files with command "'find $indir -regextype posix-extended -regex \".*(${sample_regexp}).*${READ1_REGEXP}.fastq.gz?\"'"
infiles=`find $indir -regextype posix-extended -regex ".*(${sample_regexp}).*${READ1_REGEXP}.fastq.gz?"`

echo $(date) "Going to run pipeline on " 
for f in $infiles; do 
    echo -e "\t" ${f%.fastq.gz}
done
read1=$infiles
read2=`echo $infiles | sed -e "s/${READ1_REGEXP}/${READ2_REGEXP}/g"`

##################################################
# Pipeline code
# The various steps follow a presentation on HaloPlex
# analysis obtained from Agilent
##################################################

# 1. QC
command=""
for f in $read1 $read2; do
    outdir=`dirname $f`/fastqc
    up_to_date $outdir/`basename ${f%.fastq.gz}`_fastqc/summary.txt
    if [ $? = 1 ]; then continue; fi
    mkdir -p $outdir
    cmd="$FASTQC $f -o ${outdir}"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL

# 2. Trim adapter sequence
command=""
trimfiles=""
for f in $read1; do
    trimfiles="$trimfiles ${f%.fastq.gz}.trimmed.fastq.gz"
    up_to_date ${f%.fastq.gz}.trimmed.fastq.gz
    if [ $? = 1 ]; then continue; fi
    cmd="$CUTADAPT $CUTADAPT_OPTS -a $THREEPRIME $f -o ${f%.fastq.gz}.trimmed.fastq.gz > ${f%.fastq.gz}.trimmed.fastq.cutadapt_metrics"
    command="$command\n$cmd"
done
for f in $read2; do
    trimfiles="$trimfiles ${f%.fastq.gz}.trimmed.fastq.gz"
    up_to_date ${f%.fastq.gz}.trimmed.fastq.gz
    if [ $? = 1 ]; then continue; fi
    cmd="$CUTADAPT $CUTADAPT_OPTS -a $FIVEPRIME $f -o ${f%.fastq.gz}.trimmed.fastq.gz > ${f%.fastq.gz}.trimmed.fastq.cutadapt_metrics"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL

# 3. Align sequences with bwa. Here we run command sequentially since
# bwa takes care of parallelization. From now on we run at sample level.
for f in $trimfiles; do
    echo $(date) aligning reads $f
    up_to_date ${f%.fastq.gz}.sai
    if [ $? = 1 ]; then continue; fi
    $BWA aln -t $N_CORES $BWA_REF $f > ${f%.fastq.gz}.sai 2>> $ERRFILE
done

# 4. Pair reads
sample_pfx=`for f in $read1; do echo ${f%_${READ1_REGEXP}.fastq.gz}; done`
command=""
for f in $sample_pfx; do
    label=`basename $f`
    echo $(date) pairing reads for sample $f
    up_to_date $f.sam
    if [ $? = 1 ]; then continue; fi
    echo "$BWA sampe -A -P -r \"@RG\tID:${label}\tSM:${label}\tPL:Illumina\tCN:Agilent\" $BWA_REF ${f}_${READ1_REGEXP}.trimmed.sai ${f}_${READ2_REGEXP}.trimmed.sai ${f}_${READ1_REGEXP}.trimmed.fastq.gz ${f}_${READ2_REGEXP}.trimmed.fastq.gz > $f.sam"
    cmd="$BWA sampe -A -P -r \"@RG\tID:${label}\tSM:${label}\tPL:Illumina\tCN:Agilent\" $BWA_REF ${f}_${READ1_REGEXP}.trimmed.sai ${f}_${READ2_REGEXP}.trimmed.sai ${f}_${READ1_REGEXP}.trimmed.fastq.gz ${f}_${READ2_REGEXP}.trimmed.fastq.gz > $f.sam"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL

# 5. Generate bam file
command=""
for f in $sample_pfx; do
    echo $(date) generating bam file for $f
    up_to_date $f.sort.bam
    if [ $? = 1 ]; then continue; fi
    echo "$SAMTOOLS view -bS $f.sam | $SAMTOOLS sort - $f.sort; $SAMTOOLS index $f.sort.bam "
    cmd="$SAMTOOLS view -bS $f.sam | $SAMTOOLS sort - $f.sort; $SAMTOOLS index $f.sort.bam;"
    command="$command\n$cmd"
done
echo -e $(date) $command
echo -e $command | $PARALLEL
