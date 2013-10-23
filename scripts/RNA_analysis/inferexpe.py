import os
import sys
from string import *
import math
import string
import re
import commands
import operator

if len(sys.argv) < 5:
    print """
    Usage:

    RseqQc_inferexpe.py <bed_file> <mail> <config_file> <path> 
    
    bed_file      
    mail                    eg: maya.brandi@scilifelab.se
    config_file             post_process.yaml assumes that you have specified samtools 
    version under 'custom_algorithms'/'RNA-seq analysis'
    path                    Path to analysis dir containing the tophat_out_ directories"""
    sys.exit()

bed_file        = sys.argv[1]
mail            = sys.argv[2]
config_file     = sys.argv[3]
path            = sys.argv[4]
names=commands.getoutput("ls -d tophat_out_*|sed 's/tophat_out_//g'").split('\n')
outList=[]


try:
    config  = load_config(config_file)
    extra_arg=config['sbatch']['extra_arg']
except:
    print 'ERROR: problem loading samtools version from config file'
    print  mail, bed_file, path, extra_arg

f=open("RseqQc_inferexpe.sh",'w')

print >>f, """#!/bin/bash -l
#SBATCH -A a2012043
#SBATCH -p node
#SBATCH -t 2:00:00
#SBATCH -e infer_experiment.err
#SBATCH -o infer_experiment.out
#SBATCH -J infer_experiment
#SBATCH --mail-type=ALL
#SBATCH --mail-user={0}
#SBATCH {1}""".format(mail, bed_file)


for name in names:
	DIR = str('tophat_out_'+name)
	try:
		countFile=commands.getoutput("ls "+DIR+"/"+name+".counts")
		totNum=commands.getoutput("awk '{SUM+=$2} END {print SUM}' "+countFile)
		if totNum != '':
			rRNAnum=0
			Lines=open(countFile).readlines()
			n=0
			for line in Lines:
				geneID=line.split()[0]
				if geneID in rRNAgeneList:
					n=n+1
					num=int(line.split()[1])
					rRNAnum=rRNAnum+num
			percent=round((float(rRNAnum)/int(totNum))*100,2)
			outLine=countFile.split('/')[-1].split('.')[0]+'\t'+str(percent)+'%'+'\n'
			outList.append(outLine)
        except:
                print "could not handle " + DIR
                pass	
if outList==[]:
	print 'No data found. Check count tables!'

outF=open("rRNA.quantification",'w')
outF.writelines(outList)
outF.close()
