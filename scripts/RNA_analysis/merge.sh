#!/bin/bash -l 

#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 7-00:00:00
#SBATCH -J merge
#SBATCH -e merge.err
#SBATCH -o merge.out
#SBATCH --mail-user jun.wang@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --qos=seqver

module load bioinfo-tools

module load samtools

WP=/pica/h1/funk_001/opt/scilifelab/scripts/RNA_analysis
while getopts ":p:g:m:c:d:" option; do
        case ${option} in
        p) path=${OPTARG};;
        g) gtf_file=${OPTARG};;
        m) mail=${OPTARG};;
		c) config_file=${OPTARG};;
        d) date=${OPTARG};;
        esac
done
shift $(( OPTIND - 1 ))
flowcells=''
for dir in "$@"; do
        flowcells=$flowcells" "$dir
done
flowcells=($flowcells)
nr_fc=${#flowcells[*]}
name_list=`for dir in ${flowcells[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq`
echo ${flowcells[*]}

if [ $nr_fc  -lt 2 -o $# -lt 1 ]; then
echo "
Usage:
	merge.sh -p <path> -g [gtf file] -m [mail] -c [config file]  <run dir 1> <run dir 2> ... <run dir N>

Arguments:
	<path> 
	- Path to the intermediate directory containing the runs to be analysis_$date
       	- eg: /proj/a2010002/projects/m_muurinen_11_01a/intermediate, containing 
	the three runs: 20120127B_hiseq2000 20120228A_hiseq2000 20120328A_hiseq2000

	[mail] 
	- if you want to run cufflincs and HTseq on analysis_$date samples give your mail for SLURM
	messages

	[config_file]
	- if you want to run cufflincs and HTseq on analysis_$date samples give config_file
	specifying HTseq version

       	[gtf fie]
	- Optional! If not given, the script will only merge, but not run cufflinks 
	and HTseq on the analysis_$date data.
        - Reference annotation in gtf format, used by cufflinks and HTseq

	<run dir i>
	name of directories to merge, eg
	120727_AD14B4ACXX 120727_BD1591ACXX 120823_BD16GRACXX
	

Output:
	- A new directory for the analysis_$date data, called analysis_$date and placed in the intermediate directory.
	"
else
	# merge first two flowcells
    mkdir $path/analysis_$date
    for name in ${name_list[*]};do
        samp_dir=tophat_out_$name 
        bam_file=accepted_hits_${name}.bam
        if [[ -e ${flowcells[0]}/$samp_dir && -e ${flowcells[1]}/$samp_dir ]];then ## merge read 
            echo samtools merge $path/analysis_$date/$samp_dir/$bam_file $path/${flowcells[0]}/$samp_dir/$bam_file $path/${flowcells[1]}/$samp_dir/$bam_file
            mkdir $path/analysis_$date/$samp_dir $path/analysis_$date/$samp_dir/logs
            samtools merge $path/analysis_$date/$samp_dir/$bam_file $path/${flowcells[0]}/$samp_dir/$bam_file $path/${flowcells[1]}/$samp_dir/$bam_file
            ## sum read counts
            counts_0=(`grep 'reads have been filtered out' $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
            counts_1=(`grep 'reads have been filtered out' $path/${flowcells[1]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
            counts=$((${counts_0[1]}+${counts_1[1]}))
            sorted=$((${counts_0[0]}+${counts_1[0]}))
            sed -e "s/\b${counts_0[1]}\b/$counts/g" $path/${flowcells[0]}/$samp_dir/logs/prep_reads.log|sed -e "s/\b${counts_0[0]}\b/$sorted/g" > $path/analysis_$date/$samp_dir/logs/prep_reads.log
        else if [ -e ${flowcells[0]}/tophat_out_$name ];then ## copy reads and counts
            echo "mv ${flowcells[0]}/tophat_out_$name to analysis_$date"
            mv $path/${flowcells[0]}/$samp_dir $path/analysis_$date
        else if [ -e ${flowcells[1]}/tophat_out_$name ];then ## copy reads and counts 
            echo "mv ${flowcells[1]}/tophat_out_$name to analysis_$date"
            mv $path/${flowcells[1]}/$samp_dir $path/analysis_$date
            fi;fi;fi 
    done
    # merge all remaining flowcells with "analysis_$date"
    if [ $nr_fc -gt 2 ]; then
    echo "merge all remaining flowcells with analysis_$date"
    for ((i=2; i<$nr_fc; i++));do
        for name in ${name_list[*]};do       
            samp_dir=tophat_out_$name
            bam_file=accepted_hits_${name}.bam
            if [[ -e ${flowcells[i]}/$samp_dir && -e analysis_$date/$samp_dir ]];then ## merge read
                echo samtools merge $path/analysis_$date/$samp_dir/temp_${bam_file} $path/analysis_$date/$samp_dir/$bam_file $path/${flowcells[$i]}/$samp_dir/$bam_file
				samtools merge $path/analysis_$date/$samp_dir/temp_${bam_file} $path/analysis_$date/$samp_dir/$bam_file $path/${flowcells[$i]}/$samp_dir/$bam_file
                mv $path/analysis_$date/$samp_dir/temp_${bam_file} $path/analysis_$date/$samp_dir/$bam_file
                ## sum read counts
                counts_0=(`grep 'reads have been filtered out' $path/analysis_$date/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                counts_1=(`grep 'reads have been filtered out' $path/${flowcells[$i]}/$samp_dir/logs/prep_reads.log|cut -f 1,4 -d ' '`)
                counts=$((${counts_0[1]}+${counts_1[1]}))
                sorted=$((${counts_0[0]}+${counts_1[0]}))
                sed -e "s/\b${counts_0[1]}\b/$counts/g" $path/analysis_$date/$samp_dir/logs/prep_reads.log|sed -e "s/\b${counts_0[0]}\b/$sorted/g" > $path/analysis_$date/$samp_dir/logs/temp_prep_reads.log
                mv $path/analysis_$date/$samp_dir/logs/temp_prep_reads.log $path/analysis_$date/$samp_dir/logs/prep_reads.log
            else if [ -e ${flowcells[i]}/$samp_dir ];then ## copy reads and counts
                echo "mv ${flowcells[$i]}/$samp_dir/$bam_file to analysis_$date/"
                mv $path/${flowcells[i]}/$samp_dir $path/analysis_$date
                fi;fi
            done
        done	
    fi
    rerun=`for dir in ${flowcells[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq -d`
    echo $rerun
    for i in ${rerun[*]};do
    rm $path/analysis_$date/tophat_out_${i}/accepted_hits_sorted_dupRemoved_${i}.bam
    done
    if [ "$gtf_file" != "" ]; then ## run HTseq and cufflinks on meged samples
    for i in ${rerun[*]};do
        echo make_MarkDup_HT_cuff.py $i $gtf_file $mail $path/analysis_$date $config_file
        make_MarkDup_HT_cuff.py $i $gtf_file $mail $path/analysis_$date $config_file
    done
    fi
fi
