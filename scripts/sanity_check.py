#!/usr/bin/env python

import sys
import os
import gzip
import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from operator import itemgetter, attrgetter


parser = argparse.ArgumentParser(description="Check if there are any significant contributions from unexpected organisms. " \
    "The script screens the subset of reads against the Genebank database by using BLAST version 2.2.27+ (blast/2.2.27+). " \
    "Either a gzip compressed or uncompressed fastq file is requierd. " \
    "Results are placed in the subdirectory (sanity_check) under the current folder. ")

parser.add_argument('-s', '--subset', dest="subset", default="5000", help="Number of reads for the subset (default=5000)")
parser.add_argument('-d', '--discard', dest="startFq", default="25000", help="Number of reads discard from the top (default=25000)")
parser.add_argument('-t', '--num_threads',  dest="threads", default="1", help="Number of threads (CPUs) to use in the BLAST search (default=1)")
parser.add_argument('-n', '--nhit',  dest="nhit", default="50", help="Number of top hits to report (default=50)")
parser.add_argument('file_name', help="Fastq file name, e.g. P286_109B_index2_CGATGT_L001_R1_001.fastq.gz or P286_109B_index2_CGATGT_L001_R1_001.fastq")
args = parser.parse_args()

inFile=args.file_name
subset = int(args.subset)
startFq = int(args.startFq)
lastFq = subset+startFq
threads = int(opts.threads)
nhit = int(opts.nhit)

outdir="sanity_check"
sample=inFile.split('.')[0]
tmpdir=os.path.join(outdir,"tmp")
xmldir=os.path.join(outdir,"blast_out")

Fa=sample+".sub"+str(subset)+".fa"
xml=sample+".sub"+str(subset)+".blast_out.xml"
result=sample+".sub"+str(subset)+"_sanityCheck.txt"

outFa= os.path.join(tmpdir,Fa)
outXML= os.path.join(xmldir,xml)
resultF= os.path.join(outdir,result)

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(tmpdir):
    os.mkdir(tmpdir)
if not os.path.exists(xmldir):
    os.mkdir(xmldir)

if inFile.endswith('.gz'):
    readFq=gzip.open(inFile)
else:
    readFq=open(inFile)

sub=open(outFa, 'w')
rec_no = 0
for rr in range(startFq,lastFq):
    while rec_no < rr:
        for i in range(4):
		readFq.readline()
        rec_no += 1
    anno=readFq.readline().replace("@",">")
    seq=readFq.readline()
    dump1=readFq.readline()
    dump2=readFq.readline()
    sub.write(anno)
    sub.write(seq)
    rec_no += 1
sub.close()
readFq.close()

blast_cmd=NcbiblastnCommandline(query=outFa, db="/bubo/nobackup/uppnex/blast_databases/nt", evalue=0.1, outfmt=5, out=outXML, best_hit_overhang=0.1, best_hit_score_edge=0.1, num_threads=threads, max_target_seqs=1)
blast_cmd()

result_handle = open(outXML)
blast_records= NCBIXML.parse(result_handle)

alignList=[]
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        align = str(alignment.title).split("|")[-1].split()[0:2]
        organism=' '.join(align)
        alignList.append(organism)
result_handle.close()

numAlign=len(alignList)
unique=set(alignList)
outList=[]
header=[(sample, str(startFq)+" reads discarded from the top of fastq file", str(subset)+" reads processed", str(numAlign)+" hits found", str(nhit)+" top hits reported"),("Organism", "Percentage %", "No. of hits" )]
for item in unique:
	percent=round(((float(alignList.count(item))/float(subset))*100),4)
	outList.append((item, percent, str(alignList.count(item))))
outLines=header+sorted(outList, key=itemgetter(1), reverse=True)
outF=open(resultF, 'w')
for i in range(0,(nhit+2)):
	line=list(outLines[i])
	line[1]=str(line[1])
	outF.write("%s\n" % '\t'.join(line))
outF.close()
