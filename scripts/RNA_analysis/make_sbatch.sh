#!/bin/bash -l
#SBATCH -A a2012043 
#SBATCH -p core
#SBATCH -t 01:00:00
#SBATCH -J make_sbatch
#SBATCH -e make_sbatch.err
#SBATCH -o make_sbatch.out
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
name_list=(`echo $names | tr "," "\n"`)

cd $path
for i in ${name_list[*]};do
	## read distribution
	make_RseqQc_rd.py $i $bedfile $mail $config_file $path
	sbatch RSeQC_${i}_rd.sh 

	## statistics
	get_stat.py ${i} $mail $config_file 
	sbatch ${i}_get_stat.sh
done

## quantify_rRNA
make_sbatch.py a2012043 node 01:00:00 quantify_rRNA $mail $config_file
echo "cd $path
quantify_rRNA.py $gtf_file" >> quantify_rRNA.sh
sbatch quantify_rRNA.sh

## correl
make_sbatch.py a2012043 node 01:00:00 correl $mail $config_file
echo "cd $path
R CMD BATCH '--args ${name_list[*]}' $WP/correl.R" >> correl.sh
sbatch correl.sh

## analysis report
echo "#!/bin/bash -l" > analysis_report.sh
echo "cd $path" >> analysis_report.sh
echo "analysis_report.py $run_name $project_id $names -c $config_file -r -s -d -f" >> analysis_report.sh
chmod 777 analysis_report.sh

mkdir sbatch_scripts
mv *.sh sbatch_scripts

