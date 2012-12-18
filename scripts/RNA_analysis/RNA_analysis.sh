#!/bin/bash
#	This script generates some basic statistics and a report for a RNA-seq project
#	and should be run from the intermediate directory of the project to be analysed.

if [ $# -lt 11 ]; then
  echo "Usage:
	stand in 'intermediate' and run

	RNA_analysis.sh -p <project id> -b <bed file> -g <gtf file> -m <mail> -c <config_file> <run dir 1> <run dir 2> ... <run dir N>

Arguments:
        <run dir i>
                - The name of the directory with the tophat_out_* -dirs.
		This is typically the same as the run name, such as
		20120323A_hiseq2000, but can be any name. The name of
		the run dir will also be the name set as the 'run name' 
		in the report. You might want to change this in the 
		rst file if your run dir name doesn't have an appropriate
		'run name'.

                - If more than one directory is given, the script will merge 
		the bamfiles from the diferent directories and do the analysis 
		on the merged runs

        <project id>
                - eg: M.Muurinen_11_01a

        <bed file>
                - reference gene model in bed format. Used by Ever-Seq
                to get gene body coverage and read distribution.

        <gtf fie>
                - reference annotation in gtf format, used by cufflinks and HTseq

	<mail>
		- mail adress for SLURM messages

	<config_file>
		- post_process.yaml"
  exit
fi

WP=~/opt/scilifelab/scripts/RNA_analysis
path=`pwd`
while getopts ":p:b:g:m:c:" option; do
        case ${option} in
                p) project_id=${OPTARG};;
                b) bedfile=${OPTARG};;
                g) gtf_file=${OPTARG};;
		m) mail=${OPTARG};;
		c) config_file=${OPTARG};;
        esac
done
shift $(( OPTIND - 1 ))

run_dirs=''
for dir in "$@"; do
        run_dirs=$run_dirs" "$dir
done
run_dirs=($run_dirs)
## get sample names
name_list=`for dir in ${run_dirs[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq`
names=`echo $name_list|sed -e 's/ /,/g'`
DEPENDENCY_MERGE='afterok'
DEPENDENCY_HT='afterok'
DEPENDENCY='afterok'

if [ ${#run_dirs[*]} -gt 1 ];then
	run_dir='merged'
	## megre old and new samples
	JOBID=`sbatch $WP/merge.sh -p $path ${run_dirs[*]}| sed -re 's/.+\s+([0-9]+)/\1/'`
	DEPENDENCY_MERGE=$DEPENDENCY_MERGE:$JOBID
 	echo $DEPENDENCY_MERGE
	## get names of samples to be merged
	rerun=`for dir in ${run_dirs[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq -d`
	## run HTseq and cufflinks on meged samples
	path=$path/merged
	for i in ${rerun[*]};do
		make_HT_cuff.py $i $gtf_file $mail $path $config_file	
		JOBID=`sbatch --dependency=$DEPENDENCY_MERGE "HT_cuff_$i.sh"| sed -re 's/.+\s+([0-9]+)/\1/'`
		DEPENDENCY_HT=$DEPENDENCY_HT:$JOBID
	done
	echo $DEPENDENCY_HT
	if [ $DEPENDENCY_HT = 'afterok' ]; then
		dep=" --dependency=$DEPENDENCY_MERGE"
	else 
		dep=" --dependency=$DEPENDENCY_HT"
	fi
else
	dep=""
	run_dir=${run_dirs[*]}
	path=$path/$run_dir
fi

sbatch$dep $WP/make_sbatch.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP $mail

