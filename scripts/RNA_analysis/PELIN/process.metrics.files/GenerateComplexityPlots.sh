#! /bin/bash -l

#by Pelin Akan pelin.akan@scilifelab.se
#Arguments GenerateComplexityPlots.sh

#Outputs a scatter plot of complexity curves of all samples in the project
#Requires ccurveplots.R 

space=" "
for files in `ls tophat_out_*/*.ccurve.txt | sort`; do filelist+=($files$space); done;

asize=${#filelist[@]}
echo $asize

Rscript /bubo/home/h24/mayabr/opt/scilifelab/scripts/RNA_analysis/PELIN/process.metrics.files/ccurveplots.R ${filelist[@]} 
 







