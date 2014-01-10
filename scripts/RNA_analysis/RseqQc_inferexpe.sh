#!/bin/bash 
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 2:00:00
#SBATCH -e infer_experiment.err
#SBATCH -o infer_experiment.out
#SBATCH -J infer_experiment

BED=$1

dict='{'
sep=''
for i in tophat_out_*; 
    do infer_experiment.py -i ${i}/accepted_hits_sorted_dupRemoved_*bam -r $BED  &> ${i}/RSeqC.out;
    num='"'`grep '1+-,1-+,2++,2--' ${i}/RSeqC.out|cut -f 2 -d ':'`'"';
    name='"'`echo $i|sed 's/tophat_out_//g'`'"'
    dict=${dict}${sep}${name}':'${num}
    sep=','
done
dict=${dict}'}'
echo $dict>infer_experiment.json

