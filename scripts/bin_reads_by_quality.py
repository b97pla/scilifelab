
import os
import sys
import miseq_data

PHRED_OFFSET = 64

def main():
    args = sys.argv[1:]
    
    if len(args) < 2:
        return 1
    
    process_fastq(args[0],args[1])

def print_average_quals(qualities):
    
    first = True
    for bin, avg_quality in sorted(qualities.items()):
        if first:
            print ",".join([str(i) for i in range(len(avg_quality)+1)])
            first = False
        avg_quality.insert(0,bin)
        print ",".join([str(i) for i in avg_quality])
        
def process_fastq(fastq_r1, fastq_r2, bins=[25,30], phred_offset=64):
    
    fh_r1 = miseq_data.FastQParser(fastq_r1)
    fh_r2 = miseq_data.FastQParser(fastq_r2)
    oh1 = {}
    oh2 = {}
    root1, ext1 = os.path.splitext(fastq_r1)
    root2, ext2 = os.path.splitext(fastq_r2)
    for b in bins:
        oh1[b] = miseq_data.FastQWriter("%s.Q%d%s" % (root1,b,ext1))
        oh2[b] = miseq_data.FastQWriter("%s.Q%d%s" % (root2,b,ext2))
    
    for r1 in fh_r1:
        r2 = fh_r2.next()
        assert r2[0][0:-1] == r1[0][0:-1], "FATAL: Read identifiers differ for paired reads (%s and %s)" % (r1[0],r2[0])
            
        bin = min(int(round(quality_average(r1[3].strip(),phred_offset))),int(round(quality_average(r2[3].strip(),phred_offset))))
        
        for b in bins:
            if bin >= b:
                oh1[b].write(r1)
                oh2[b].write(r2)
        
    for oh in oh1.values() + oh2.values():
        oh.close()


def quality_average(quality_str, phred_offset):
    sum = 0.
    for qual in quality_str:
        sum += ord(qual)
    sum -= len(quality_str)*phred_offset
    return sum/len(quality_str)
    
if __name__ == "__main__":
    sys.exit(main())
    