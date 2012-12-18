import sys
from bcbio.pipeline.config_loader import load_config

if len(sys.argv) < 6:
        print """
Usage:

make_RseqQc_gbc.py  <sample name> <bed_file> <mail> <config_file> <path>

        <sample name>           This name: /tophat_out_<sample name>
        <bed_file>      
        <mail>                  eg: maya.brandi@scilifelab.se
        <config_file>           post_process.yaml assumes that you have specified samtools 
                                version under 'custom_algorithms'/'RNA-seq analysis'
        <path>                  Path to analysis dir containing the tophat_out_ directories
        """
        sys.exit()

name            = sys.argv[1]
bed_file        = sys.argv[2]
mail            = sys.argv[3]
config_file     = sys.argv[4]
path            = sys.argv[5]
try:
        config  = load_config(config_file)
        tools   = config['custom_algorithms']['RNA-seq analysis']
        sam     = tools['sam']+'/'+tools['sam_version']
except:
        print 'ERROR: problem loading samtools version from config file'


f=open("RSeQC_"+name+"_gbc.sh",'w')

print >>f, """#!/bin/bash -l
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 50:00:00
#SBATCH --qos=seqver
#SBATCH -e RSeQC_gbc_{0}.err
#SBATCH -o RSeQC_gbc_{0}.out
#SBATCH -J RSeQC_gbc_{0}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={1}

module unload samtools
module load {3}

cd {4}

samtools view /tophat_out_{0}/accepted_hits_sorted_dupRemoved_{0}.bam | geneBody_coverage.py -i - -r {2} -o {0}
CMD BATCH {0}.geneBodyCoverage_plot.r
mv geneBody_coverage.pdf {0}_geneBody_coverage.pdf""".format(name, mail, bed_file, sam, path)
