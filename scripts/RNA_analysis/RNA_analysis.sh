#!/bin/bash -l
#	This script generates some basic statistics and a report for a RNA-seq project
#	and should be run from the intermediate directory of the project to be analysed.

while getopts ":p:b:g:m:c:e:a:s:d:f:" option; do
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
        f) single=${OPTARG};;
    esac
done
shift $(( OPTIND - 1 ))

run_dirs=''
for dir in "$@"; do
    run_dirs=$run_dirs" "$dir
done

WP=/pica/h1/funk_001/opt/scilifelab/scripts/RNA_analysis
path=`pwd`
run_dir=analysis_$date
analysis_path=$path/$run_dir
mkdir $analysis_path
run_dirs=($run_dirs)
name_list=`for dir in ${run_dirs[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq`
for dir in ${run_dirs[*]};do echo $dir >> $analysis_path/analysis_log.txt;ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g' >> $analysis_path/analysis_log.txt;done
names=`echo $name_list|sed -e 's/ /,/g'`

if [ ${#run_dirs[*]} -gt 1 ];then
	## megre old and new samples
	JOBID=`sbatch $WP/merge.sh -p $path -d $date ${run_dirs[*]}| sed -re 's/.+\s+([0-9]+)/\1/'`
	merge_dep=" --dependency=afterok:$JOBID"
	## get names of samples that are being merged
	rerun=`for dir in ${run_dirs[*]};do ls -d $dir/tophat_out_*|cut -f 2 -d '/'|sed 's/tophat_out_//g';done|sort|uniq -d`
else
    mv $path/${run_dirs[0]}/tophat_out* $analysis_path
    merge_dep=""
fi

if [ $analysis = 'False' ]; then
    analysisi_list=${rerun[*]}
else
    analysisi_list=$name_list
fi

if [ $single = 'False' ]; then
    single=''
else
    single='-e'
fi

DEPENDENCY='afterok'
for i in $analysisi_list;do
    make_MarkDup_HT_cuff.py $i $gtf_file $mail $analysis_path $config_file $stranded $extra_arg
    JOBID=`sbatch $merge_dep "MarkDup_HT_cuff_$i.sh"| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEPENDENCY=$DEPENDENCY:$JOBID
done
if [ $DEPENDENCY = 'afterok' ]; then
    dep=""
else 
    dep=" --dependency=$DEPENDENCY"
    echo "If you get a warning telling that the dependensy list is to long, you will have to restart the seting_dependensies.sh script manualy. This should be done after the folowing jobs are finished:
    
    $DEPENDENCY
    
    When all the jobbs are finished, start the rest of the pipeline by giving the command:

    sbatch $WP/seting_dependensies.sh $names $bedfile $project_id $config_file $run_dir $analysis_path $gtf_file $WP $mail $single

    This command is saved in the text file seting_dependensies.txt 
    "
    echo sbatch $WP/seting_dependensies.sh $names $bedfile $project_id $config_file $run_dir $analysis_path $gtf_file $WP $mail $single  >> seting_dependensies.txt
fi

sbatch$dep $WP/seting_dependensies.sh $names $bedfile $project_id $config_file $run_dir $analysis_path $gtf_file $WP $mail $single
