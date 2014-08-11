#! /bin/bash -l
#SBATCH -A a2012043 
#SBATCH -p core
#SBATCH -t 01:00:00
#SBATCH -J complexity
#SBATCH -e complexity.err
#SBATCH -o complexity.out

#by Francesco Vezzi francesco.vezzi@scilifelab.se
#Arguments GenerateComplexityPlots.sh

#Outputs a scatter plot of complexity curves of all samples in the project
#Requires plot_complexity_curves.py 
wp=$1

unset num
unset filelist
##List all *.ccurve.txt files, and get the number of the 'total_reads' in *.ccurve.txt for all samples
space=" "
for files in `ls tophat_out_*/*.ccurve.txt | sort`; do filelist+=($files$space); n=`tail -n 1 ${files}|cut -f 1`; num+=($n$space); done;
asize=${#filelist[@]}
##set x-max value, should be min(total_reads) of all samples
FLOAT=`echo ${num[@]} |sed 's/ /\n/g' |sort -n|sed -n '1p'`
xmax=${FLOAT/.*}
python ${wp}/plot_complexity_curves.py --ccurves ${filelist[@]} --x-max ${xmax} 
