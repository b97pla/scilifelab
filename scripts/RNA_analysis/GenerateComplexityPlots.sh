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
space=" "
for files in `ls tophat_out_*/*.ccurve.txt | sort`; do filelist+=($files$space); n=`tail -n 1 ${files}|cut -f 1`; num+=($n$space); done;
asize=${#filelist[@]}
xmax=`echo ${num[@]} |sed 's/ /\n/g' |sort -n|sed -n '1p'`

python ${wp}/plot_complexity_curves.py --ccurves ${filelist[@]} --x-max ${xmax} 
