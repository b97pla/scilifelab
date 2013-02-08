"""
Reads a FastQ file from stdin and writes a file with unique records to sdout
usage:
    %s < in.fastq > out.unique.fastq

see: http://hackmap.blogspot.com/2010/10/bloom-filter-ing-repeated-reads.html
"""
from bloomfaster import Elf
import collections
import sys

# Slight modification to read from input file instead of stdin
from scilifelab.utils.fastq_utils import (FastQParser, FastQWriter)

__doc__ %= sys.argv[0]
if len(sys.argv) > 2:
    print sys.argv
    print __doc__
    sys.exit()

print >>sys.stderr, "Command: ", " ".join(sys.argv)
infile = sys.argv[1]
fp = FastQParser(infile)
for _ in fp:
    pass
records = fp.rread()
print >>sys.stderr, records, "records in file ", infile

# say 1 out of 1000 is false positive.
bloom = Elf(records, error_rate=1e-3)
fp.seek(0)
checks = []
for _,seq,_,_ in fp:
    if seq in bloom:
        checks.append(seq)
    bloom.add(seq)

# now checks contains anything that could be a duplicate according to
# the bloomfilter. for some, they were false positives.
# for actual duplicated, just choose the first, but can also sort by quality.
fp.seek(0)
checks = frozenset(checks)
print >>sys.stderr, "checking %s potential duplicates in a python set" \
                                            % len(checks)
outfile = "%s-unique.fastq.gz" % infile.split(".")[0]                                            
fw = FastQWriter(outfile)

putative_false_positives = collections.defaultdict(int)
for header,seq,plus,qual in fp:
    # it appeared only once, so just print it.
    if not seq in checks:
        fw.write([header,seq,plus,qual])
        continue
    # it appeared in the bloom-filter > 1x, so track and make sure not
    # to print any others with same sequence.
    putative_false_positives[seq] += 1
    if putative_false_positives[seq] > 1:
        continue
    fw.write([header,seq,plus,qual])

fw.close()
false_positives = sum(1 for count in putative_false_positives.values() \
                                                        if count == 1)
print >>sys.stderr,  false_positives, "false-positive duplicates in the bloom filter"
