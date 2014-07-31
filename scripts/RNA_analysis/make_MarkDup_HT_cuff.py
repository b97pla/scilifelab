#!/usr/bin/env python

import sys
from bcbio.pipeline.config_loader import load_config
if len(sys.argv) < 6:
    print """Usage:

make_MarkDup_HT_cuff.py	 <sample_name> <gtf_file> <mail> <path> <config_file> <stranded> [extra_arg]

        sample name		This name: /tophat_out_<sample name>
	gtf_file	
	mail			eg: jun.wang@scilifelab.se
	path			path to analysis directory containing the tophat_out_dirs
	config_file		post_process.yaml assumes that you have specfied cufflinks 
				and HT-seq versions under 'custom_algorithms'/'RNA-seq analysis'
    stranded        True/False
        extra_arg		extra SBATCH argumet for high priority job
        """
    sys.exit()


name	 	= sys.argv[1]
gtf_file 	= sys.argv[2]
mail	 	= sys.argv[3]
path     	= sys.argv[4]
config_file	= sys.argv[5]
stranded    = sys.argv[6]
try:
    extra_arg = "#SBATCH "+sys.argv[7]
except:
    extra_arg = ""
    pass

try:
    config = load_config(config_file)
    extra_arg = config['sbatch']['extra_arg']
    tools = config['custom_algorithms']['RNA-seq analysis']
    picard_tools = tools['picard_tools']
    picard_version = tools['picard_version']
    quant = tools['quantifyer']+'/'+tools['quantifyer_version']
    counts = tools['counts']+'/'+tools['counts_version']
    rseqc_version = tools['rseqc_version']
    if stranded == 'True':
        aligner_libtype = tools['aligner_libtype']
        ht_stranded = 'reverse'
    else:
        aligner_libtype = ''
        ht_stranded = 'no'
except: 
    print 'ERROR: problem loading cufflinks and HT-seq versions from config file'

tophat_out_path = "{0}/tophat_out_{1}".format(path,name) 
f=open("MarkDup_HT_cuff_"+name+".sh",'w')
print >>f, """#!/bin/bash -l 
#SBATCH -A a2012043
#SBATCH -p core -n 8
#SBATCH -t 40:00:00
#SBATCH -e MarkDup_HT_cuff_{0}.err
#SBATCH -o MarkDup_HT_cuff_{0}.out
#SBATCH -J MarkDup_HT_cuff_{0}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={3}
#SBATCH {7}

module load bioinfo-tools
module load {4}
module load {5}
module load picard/{9}

cd {1}
java -Xmx2g -jar {8}/SortSam.jar INPUT=accepted_hits_{0}.bam OUTPUT=accepted_hits_sorted_{0}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -jar {8}/MarkDuplicates.jar INPUT=accepted_hits_sorted_{0}.bam OUTPUT=accepted_hits_sorted_dupRemoved_{0}.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE={0}_picardDup_metrics VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -jar {8}/CollectInsertSizeMetrics.jar INPUT=accepted_hits_{0}.bam OUTPUT={0}.picard_estimated_insert_size HISTOGRAM_FILE={0}.picard_histogram_file

samtools view accepted_hits_sorted_dupRemoved_{0}.bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam
python -m HTSeq.scripts.count -s {12} -q accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam {2} > {0}.counts_dup_rem
rm accepted_hits_sorted_dupRemoved_prehtseq_{0}.sam

samtools view accepted_hits_sorted_{0}.bam |sort > accepted_hits_sorted_prehtseq_{0}.sam
python -m HTSeq.scripts.count -s {12} -q accepted_hits_sorted_prehtseq_{0}.sam {2} > {0}.counts
rm accepted_hits_sorted_prehtseq_{0}.sam

samtools index accepted_hits_sorted_{0}.bam
cufflinks -p 8 {11} -G {2} -o cufflinks_out_{0} accepted_hits_sorted_{0}.bam
""".format(name, tophat_out_path, gtf_file, mail, quant, counts, path, extra_arg, picard_tools,picard_version,rseqc_version,aligner_libtype,ht_stranded)


