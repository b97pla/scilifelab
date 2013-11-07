#!/bin/bash -l
#	This script generates some basic statistics and a report for a RNA-seq project
#	and should be run from the intermediate directory of the project to be analysed.

WP=~/opt/scilifelab/scripts/RNA_analysis
path=`pwd`

while getopts ":p:b:g:m:c:e:a:s:d:" option; do
    case ${option} in
        p) project_id=${OPTARG};;   
        b) bedfile=${OPTARG};;
        g) gtf_file=${OPTARG};;
        m) mail=${OPTARG};;
        c) config_file=${OPTARG};;
        e) extra_arg=${OPTARG};;
        a) analysis=${OPTARG};;
        s) stranded=${OPTARG};;
        d) date=${OPTARG};;
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
if [ ${#run_dirs[*]} -gt 1 ];then
	run_dir='merged_'$date
	## megre old and new samples
	JOBID=`sbatch $WP/merge.sh -p $path -d $date ${run_dirs[*]}| sed -re 's/.+\s+([0-9]+)/\1/'`
	merge_dep=" --dependency=afterok:$JOBID"
	## get names of samples that are being merged
	rerun=`for dir in ${run_dirs[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq -d`
	path=$path/merged_$date
else
    merge_dep=""
    run_dir=${run_dirs[*]}
    path=$path/$run_dir
fi

if [ $analysis = 'False' ]; then
    analysisi_list=${rerun[*]}
else
    analysisi_list=$name_list
fi

DEPENDENCY='afterok'
for i in $analysisi_list;do
    make_MarkDup_HT_cuff.py $i $gtf_file $mail $path $config_file $stranded $extra_arg
    JOBID=`sbatch $merge_dep "MarkDup_HT_cuff_$i.sh"| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEPENDENCY=$DEPENDENCY:$JOBID
done
if [ $DEPENDENCY = 'afterok' ]; then
    dep=""
else 
    dep=" --dependency=$DEPENDENCY"
fi
echo $dep $WP/make_sbatch.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP $mail
#sbatch$dep $WP/make_sbatch.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP $mail
sbatch$dep $WP/seting_dependensies.sh $names $bedfile $project_id $config_file $run_dir $path $gtf_file $WP $mail
