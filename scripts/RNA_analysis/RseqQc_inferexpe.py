import sys
from bcbio.pipeline.config_loader import load_config

if len(sys.argv) < 5:
    print """
Usage:

RseqQc_inferexpe.py <bed_file> <mail> <config_file> <path> 

        bed_file      
        mail                  	eg: maya.brandi@scilifelab.se
        config_file           	post_process.yaml assumes that you have specified samtools 
                                version under 'custom_algorithms'/'RNA-seq analysis'
        path                  	Path to analysis dir containing the tophat_out_ directories"""
    sys.exit()

bed_file        = sys.argv[1]
mail            = sys.argv[2]
config_file     = sys.argv[3]
path            = sys.argv[4]

try:
    config  = load_config(config_file)
    extra_arg = config['sbatch']['extra_arg']
    rseqc_version = config['custom_algorithms']['RNA-seq analysis']['rseqc_version']
except:
    print 'ERROR: problem loading samtools version from config file'

f=open("RseqQc_inferexpe.sh",'w')

print >>f, """#!/bin/bash -l 
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 2:00:00
#SBATCH -e infer_experiment.err
#SBATCH -o infer_experiment.out
#SBATCH -J infer_experiment
#SBATCH --mail-type=ALL
#SBATCH --mail-user="""+mail+"""
#SBATCH """+extra_arg+"""

module load bioinfo-tools
module load rseqc/"""+rseqc_version+"""

cd """+path+"""

dict='{'
sep=''
for i in tophat_out_*;do 
    infer_experiment.py -i ${i}/accepted_hits_*bam -r """+bed_file+"""  &> ${i}/RSeqC.out;
    num='"'`grep '1+-,1-+,2++,2--' ${i}/RSeqC.out|cut -f 2 -d ':'`'"';
    name='"'`echo ${i}|sed 's/tophat_out_//g'`'"'
    dict=${dict}${sep}${name}':'${num}
    sep=','
done
dict=${dict}'}'
echo ${dict}>infer_experiment.json"""


