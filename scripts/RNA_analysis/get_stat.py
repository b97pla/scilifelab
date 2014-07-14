import sys
from bcbio.pipeline.config_loader import load_config


if len(sys.argv) < 3:        
    print """
Usage:

get_stat.py  <sample_name> <mail> [config]
    sample_name		This name: /tophat_out_<sample name>
    mail                  	eg: jun.wang@scilifelab.se
    config            	config file"""
    sys.exit()

name	  = sys.argv[1]
mail	  = sys.argv[2]


try:
    config    = load_config(sys.argv[3])
    extra_arg = "#SBATCH "+config['sbatch']['extra_arg']
    rseqc_version = config['custom_algorithms']['RNA-seq analysis']['rseqc_version']
except:
    extra_arg = ""
    rseqc_version = ""
    pass

f=open(name+"_get_stat.sh",'w')
print >>f, """#! /bin/bash -l

#SBATCH -A a2012043
#SBATCH -p core
#SBATCH -t 2:00:00
#SBATCH -J get_stat"""+ name + """
#SBATCH -o get_stat"""+ name + """.out
#SBATCH -e get_stat"""+ name + """.err
#SBATCH --mail-user """+ mail +"""
#SBATCH --mail-type=ALL
""" + extra_arg + """

module add bioinfo-tools
module load samtools
module load rseqc/""" + rseqc_version + """

NAME="""+name+"""

BAM=tophat_out_${NAME}/accepted_hits_${NAME}.bam
DUP_REM_BAM=tophat_out_${NAME}/accepted_hits_sorted_dupRemoved_${NAME}.bam
bam_stat.py --version >> get_stat${NAME}.err
echo '' >> get_stat${NAME}.err
echo '===== Mapping statistics based on all reads =====' >> get_stat${NAME}.err
bam_stat.py -i $BAM
echo '' >> get_stat${NAME}.err
echo '===== Mapping statistics based on reads after duplicate removal =====' >> get_stat${NAME}.err
bam_stat.py -i $DUP_REM_BAM
SPLICED=`samtools view $DUP_REM_BAM |cut -f6 |grep -c N`
echo '===== Number of spliced alignments (dupRem) ====='>> get_stat${NAME}.err
echo 'spliced:          '$SPLICED >> get_stat${NAME}.err

mv get_stat${NAME}.err tophat_out_${NAME}/stat${NAME}

"""
