import sys
from bcbio.pipeline.config_loader import load_config
if len(sys.argv) < 6:
        print """
Usage:

make_HT_cuff.py	 <sample_name> <gtf_file> <mail> <path> <config_file> [extra_arg]

        sample name		This name: /tophat_out_<sample name>
	gtf_file	
	mail			eg: maya.brandi@scilifelab.se
	path			path to analysis directory containing the tophat_out_dirs
	config_file		post_process.yaml assumes that you have specfied cufflinks 
				and HT-seq versions under 'custom_algorithms'/'RNA-seq analysis'
        extra_arg		extra SBATCH argumet for high priority job
        """
	sys.exit()


name	 	= sys.argv[1]
gtf_file 	= sys.argv[2]
mail	 	= sys.argv[3]
path     	= sys.argv[4]
config_file	= sys.argv[5]
try:
        extra_arg = "#SBATCH "+sys.argv[6]
except:
        extra_arg = ""
        pass

try:
        config   = load_config(config_file)
	extra_arg= config['sbatch']['extra_arg']
        tools    = config['custom_algorithms']['RNA-seq analysis']
        quant    = tools['quantifyer']+'/'+tools['quantifyer_version']
        counts   = tools['counts']+'/'+tools['counts_version']
except: 
        print 'ERROR: problem loading cufflinks and HT-seq versions from config file'

tophat_out_path = "{0}/tophat_out_{1}".format(path,name) 
f=open("HT_cuff_"+name+".sh",'w')
print >>f, """#!/bin/bash -l
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 10:00:00
#SBATCH -e HT_cuff_{0}.err
#SBATCH -o HT_cuff_{0}.out
#SBATCH -J HT_cuff_{0}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={3}
#SBATCH {7}

module unload cufflinks
module load {4}
module load {5}

samtools view {1}/accepted_hits_sorted_dupRemoved_{0}.bam |sort > {1}/accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam
python -m HTSeq.scripts.count -s no -q {1}/accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam {2} > {1}/{0}.counts
rm {1}/accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam
samtools index {1}/accepted_hits_sorted_dupRemoved_{0}.bam
cufflinks -p 8 -G {2} -o {1}/cufflinks_out_{0} {6}/tophat_out_{0}/accepted_hits_sorted_dupRemoved_{0}.bam""".format(name, tophat_out_path, gtf_file, mail, quant, counts, path, extra_arg)

