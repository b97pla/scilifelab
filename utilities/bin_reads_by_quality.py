
import os
import sys
import fastq_utils
import argparse

PHRED_OFFSET = 64

def main():
    
    parser = argparse.ArgumentParser(description="Filter reads from a pair of FastQ files based on the average quality."\
                                     "If the average quality of one of the reads in the pair is below the given threshold, "\
                                     "the pair is discarded. Output is a file named as INPUT.Q[T].[EXT], where T is the threshold "\
                                     "and EXT is the file extension. Accepts uncompressed or gzip-compressed input files")

    parser.add_argument('-T','--threshold', action='store', default=20, 
                        help="if any read in the pair has an average quality below this threshold, the pair is discarded. Default is 20.")
    parser.add_argument('-p','--phred', action='store', default=33, 
                        help="the Phred quality score offset. Default is 33 (Sanger)")
    parser.add_argument('fastq1', action='store', default=None, 
                        help="the first sequence file of the pair")
    parser.add_argument('fastq2', action='store', default=None, 
                        help="the second sequence file of the pair")
    
    args = parser.parse_args()
    process_fastq(args.fastq1, args.fastq2, [int(args.threshold)], int(args.phred))

def print_average_quals(qualities):
    
    first = True
    for bin, avg_quality in sorted(qualities.items()):
        if first:
            print ",".join([str(i) for i in range(len(avg_quality)+1)])
            first = False
        avg_quality.insert(0,bin)
        print ",".join([str(i) for i in avg_quality])
        
def process_fastq(fastq_r1, fastq_r2, bins, phred_offset):
    
    fh_r1 = fastq_utils.FastQParser(fastq_r1)
    fh_r2 = fastq_utils.FastQParser(fastq_r2)
    oh1 = {}
    oh2 = {}
    root1, ext1 = os.path.splitext(fastq_r1)
    root2, ext2 = os.path.splitext(fastq_r2)
    for b in bins:
        oh1[b] = fastq_utils.FastQWriter("%s.Q%d%s" % (root1,b,ext1))
        oh2[b] = fastq_utils.FastQWriter("%s.Q%d%s" % (root2,b,ext2))
    
    for r1 in fh_r1:
        r2 = fh_r2.next()
        r1h = r1[0].split()
        r2h = r2[0].split()
        assert r2h[0] == r1h[0] and r2h[1][1:] == r1h[1][1:], "FATAL: Read identifiers differ for paired reads (%s and %s)" % (r1[0],r2[0])
            
        bin = min(int(round(fastq_utils.avgQ(r1,phred_offset))),int(round(fastq_utils.avgQ(r2,phred_offset))))
        
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
    