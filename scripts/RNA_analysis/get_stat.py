import sys

if len(sys.argv) < 3:        
	print """
Usage:

get_stat.py  <sample name> <mail>

        <sample name>           This name: /tophat_out_<sample name>
        <mail>                  eg: maya.brandi@scilifelab.se
        """
        sys.exit()

name	= sys.argv[1]
mail	= sys.argv[2]

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
#SBATCH --qos=seqver

module add bioinfo-tools
module load samtools

NAME="""+name+"""

BAM=tophat_out_${NAME}/accepted_hits_${NAME}.bam
DUP_REM_BAM=tophat_out_${NAME}/accepted_hits_sorted_dupRemoved_${NAME}.bam
echo '' >> get_stat${NAME}.err
echo '===== Mapping statistics based on all reads =====' >> get_stat${NAME}.err
bam_stat.py -i $BAM >> tophat_out_${NAME}/stat_${NAME}
echo '' >> get_stat${NAME}.err
echo '===== Mapping statistics based on reads after duplicate removal =====' >> get_stat${NAME}.err
bam_stat.py -i $DUP_REM_BAM >> tophat_out_${NAME}/stat_${NAME}
SPLICED=`samtools view $DUP_REM_BAM |cut -f6 |grep -c N`
echo '===== Number of spliced alignments (dupRem) ====='>> get_stat${NAME}.err
echo 'spliced:          '$SPLICED >> get_stat${NAME}.err

mv get_stat${NAME}.err tophat_out_${NAME}/stat${NAME}

"""
