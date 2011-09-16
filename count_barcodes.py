
# This script uses a few assumptions about how the barcoding works.
# Since the aim is to look at Illumina HiSeq barcodes, the assumption is that the barcode
# is located at the 3' end of read 1 for each paired-end read. 

import sys
from operator import itemgetter


if len(sys.argv)<3:
    print "python", sys.argv[0], "<fastq file containing barcode sequence> <position where bar code starts>"
    sys.exit(0)
    
bcodes = {} # collect counts for all observed barcodes

cntr = 0
last_was_header = False
pos = int(sys.argv[2])

for line in open(sys.argv[1]):
    if not last_was_header:
        if line[0]=="@": last_was_header=True
        continue

    bcode = line[pos:].strip()
    if bcodes.has_key(bcode): bcodes[bcode] += 1
    else: bcodes[bcode] = 1
    if last_was_header: last_was_header = False

#for k in bcodes.keys():
#    print k, bcodes[k]

for e in sorted(bcodes.items(), key=itemgetter(1)):
    print e[0] + "\t" + str(e[1])

