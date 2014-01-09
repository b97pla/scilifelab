#!/bin/bash -l
#SBATCH -A a2012043 
#SBATCH -p core
#SBATCH -t 01:00:00
#SBATCH -J seting_dependensies
#SBATCH -e seting_dependensies.err
#SBATCH -o seting_dependensies.out
while getopts ":e:" option; do
        case ${option} in
                e) extra_arg=${OPTARG};;
        esac
done
shift $(( OPTIND - 1 ))

names=$1
bedfile=$2
project_id=$3
config_file=$4
run_name=$5
path=$6
gtf_file=$7
WP=$8
mail=$9
single=${10}
name_list=(`echo $names | tr "," "\n"`)

cd $path
DEPENDENCY='afterok'
DEP_REPORT='afterok'
for i in ${name_list[*]};do
	## read distribution
    echo 'read distribution'
	make_RseqQc_rd.py $i $bedfile $mail $config_file $path
    echo "..."
	JOB=`sbatch RSeQC_${i}_rd.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEP_REPORT=$DEP_REPORT:$JOB
    echo "make_RseqQc_gbc.py" 

    make_RseqQc_gbc.py $i $bedfile $mail $config_file $path
    echo '....'
    JOB=`sbatch RSeQC_${i}_gbc.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEP_REPORT=$DEP_REPORT:$JOB
	## statistics
    echo 'statistics'
	get_stat.py ${i} $mail $config_file 
    echo '...'
	JOB=`sbatch ${i}_get_stat.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEP_REPORT=$DEP_REPORT:$JOB
    echo 'complexity'
    ## lib complexity
    make_complexity_plots.py ${i} $mail $config_file $path
    JOBID=`sbatch complexity_${i}.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
    DEPENDENCY=$DEPENDENCY:$JOBID
done
echo $DEPENDENCY
## make complexity plot
JOB=`sbatch --dependency=$DEPENDENCY ${WP}/GenerateComplexityPlots.sh ${WP}| sed -re 's/.+\s+([0-9]+)/\1/'`
DEP_REPORT=$DEP_REPORT:$JOB

## seqQc_inferexpe
RseqQc_inferexpe.py $bedfile $mail $config_file $path
JOB=`sbatch RseqQc_inferexpe.sh $bedfile| sed -re 's/.+\s+([0-9]+)/\1/'`
DEP_REPORT=$DEP_REPORT:$JOB

## quantify_rRNA
make_sbatch.py a2012043 node 01:00:00 quantify_rRNA $mail $config_file
echo "cd $path
quantify_rRNA.py $gtf_file" >> quantify_rRNA.sh
JOB=`sbatch quantify_rRNA.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
DEP_REPORT=$DEP_REPORT:$JOB

## correl
make_sbatch.py a2012043 node 01:00:00 correl $mail $config_file
echo "cd $path
R CMD BATCH '--args ${name_list[*]}' $WP/correl.R" >> correl.sh
JOB=`sbatch correl.sh| sed -re 's/.+\s+([0-9]+)/\1/'`
DEP_REPORT=$DEP_REPORT:$JOB

## analysis report
cp $WP/sll_logo.gif .
make_sbatch.py a2012043 core 01:00:00 analysis_report $mail $config_file
echo "cd $path
analysis_report.py $project_id -c $config_file -r -s -d -f -g -w -b $single" >> analysis_report.sh
sbatch --dependency=$DEP_REPORT analysis_report.sh
mkdir sbatch_scripts
mv *.sh sbatch_scripts

