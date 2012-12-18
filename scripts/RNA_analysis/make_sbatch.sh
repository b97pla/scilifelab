#!/bin/bash -l
#SBATCH -A a2012043 
#SBATCH -p core
#SBATCH -t 01:00:00
#SBATCH -J make_sbatch
#SBATCH -e make_sbatch.err
#SBATCH -o make_sbatch.out

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

## gene body coverage
#make_RseqQc_gbc.py $i $bedfile $mail $config_file $path

## statistics
get_stat.py ${i} $mail
sbatch ${i}_get_stat.sh
done


## quantify_rRNA
echo "#! /bin/bash -l" > quantify_rRNA.sh
echo "#SBATCH -A a2012043" >> quantify_rRNA.sh
echo "#SBATCH -p node" >> quantify_rRNA.sh
echo "#SBATCH -t 15:00:00" >> quantify_rRNA.sh
echo "#SBATCH -J quantify_rRNA" >> quantify_rRNA.sh
echo "#SBATCH -e quantify_rRNA.err" >> quantify_rRNA.sh
echo "#SBATCH -o quantify_rRNA.out" >> quantify_rRNA.sh
echo "#SBATCH --mail-user $mail" >> quantify_rRNA.sh
echo "#SBATCH --mail-type=ALL" >> quantify_rRNA.sh
echo "cd $path
quantify_rRNA.py $gtf_file" >> quantify_rRNA.sh

sbatch quantify_rRNA.sh


## correl
echo "#! /bin/bash -l" > correl.sh
echo "#SBATCH -A a2012043" >> correl.sh
echo "#SBATCH -p node" >> correl.sh
echo "#SBATCH -t 01:00:00" >> correl.sh
echo "#SBATCH -J correl" >> correl.sh
echo "#SBATCH -e correl.err" >> correl.sh
echo "#SBATCH -o correl.out" >> correl.sh
echo "#SBATCH --mail-user $mail" >> correl.sh
echo "#SBATCH --mail-type=ALL" >> correl.sh
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

