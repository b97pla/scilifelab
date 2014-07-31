import sys
from bcbio.pipeline.config_loader import load_config
if len(sys.argv) < 5:
    print """
Usage:

make_complexity_plots.py  <sample_name> <mail> <config_file> <path>

        sample_name		This name: /tophat_out_<sample name>
        mail                  	eg: jun.wang@scilifelab.se
        config_file           	post_process.yaml assumes that you have specified samtools
                                version under 'custom_algorithms'/'RNA-seq analysis'
        path                  	Path to analysis dir containing the tophat_out_ directories
        """
    sys.exit()

name            = sys.argv[1]
mail            = sys.argv[2]
config_file     = sys.argv[3]
path            = sys.argv[4]

try:
    config  = load_config(config_file)
    extra_arg = config['sbatch']['extra_arg']
    tools   = config['custom_algorithms']['RNA-seq analysis']
    bam     = tools['bamtools']+'/'+tools['bamtools_version']
    preseq  = tools['Preseq']+'/'+tools['Preseq_version']
except:
    print 'ERROR: problem loading samtools version from config file'


f=open("complexity_"+name+".sh",'w')

print >>f, """#!/bin/bash -l
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 10:00:00
#SBATCH -e complexity_{0}.err
#SBATCH -o complexity_{0}.out
#SBATCH -J complexity_{0}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={1}
#SBATCH {5}

module load bioinfo-tools
module load {3}
module load {2}

cd {4}

preseq lc_extrap -v -B tophat_out_{0}/accepted_hits_sorted_{0}.bam -o tophat_out_{0}/{0}.ccurve.txt

""".format(name, mail, preseq, bam, path, extra_arg)
