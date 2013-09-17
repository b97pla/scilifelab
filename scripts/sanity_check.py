#!/usr/bin/env python

import sys
import os
import gzip
import random
import optparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from operator import itemgetter, attrgetter

usage="""
    Check if there are any significant contributions from unexpected organisms.
    The script screens the subset of reads against the Genebank database by using BLAST version 2.2.27+ (blast/2.2.27+)

Usage:
        sanity_check.py <*.fastq.gz or *.fastq>
    The first option is mandatory, and there are further flages you can use:

-s, --subset: number of reads for the subset (default=5000)
-d, --discard: number of reads discard from the top (default=25000)
-t, --num_threads: number of threads (CPUs) to use in the BLAST search (default=1)
-n, --nhit: number of top hits to report (default=50)

Output:
    a new subdirectory "sanity_check" under the current folder.

"""

if len(sys.argv) < 2:
        sys.exit(usage)

inFile = sys.argv[1]

parser = optparse.OptionParser()
parser.add_option('-s', '--subset', type="int", dest="subset", default="5000", help="Number of reads for the subset (default=5000)")
parser.add_option('-d', '--discard', type="int", dest="startFq", default="25000", help="Number of reads discard from the top (default=25000)")
parser.add_option('-t', '--num_threads', type="int", dest="threads", default="1", help="Number of threads (CPUs) to use in the BLAST search (default=1)")
parser.add_option('-n', '--nhit', type="int", dest="nhit", default="50", help="Number of top hits to report (default=50)")
(opts, args) = parser.parse_args()

subset = opts.subset
startFq =opts.startFq
lastFq = subset*2+startFq
threads = opts.threads
nhit = opts.nhit

outdir="sanity_check"
sample=inFile.split('.')[0]
outFa= outdir+"/tmp/"+sample+".sub"+str(subset)+".fa"
outFq= outdir+"/tmp/"+sample+".sub"+str(subset)+".fastq"
outXML= outdir+"/blast_out/"+sample+".sub"+str(subset)+".blast_out.xml"
resultF= outdir+"/"+sample+".sub"+str(subset)+"_sanityCheck.txt"

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(outdir+"/tmp"):
    os.mkdir(outdir+"/tmp")
if not os.path.exists(outdir+"/blast_out"):
    os.mkdir(outdir+"/blast_out")

if inFile.endswith('.gz'):
    readFq=gzip.open(inFile)
    rand_records = sorted(random.sample(xrange(startFq,lastFq), subset))
else:
    readFq=open(inFile)
    rand_records = sorted(random.sample(xrange(startFq,lastFq), subset))

sub=open(outFq, 'w')
rec_no = -1
for rr in rand_records:
    while rec_no < rr:
        for i in range(4):readFq.readline()
        rec_no += 1
    for i in range(4):
        sub.write(readFq.readline())
    rec_no += 1
sub.close()
readFq.close()

read = SeqIO.parse(outFq,"fastq")
handle=open(outFa, "w")
for record in read:
	handle.write(">"+record.id+"\n"+str(record.seq)+"\n")
handle.close()

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

