
import os
import sys
import miseq_data

PHRED_OFFSET = 33

def main():
    args = sys.argv[1:]
    
    if len(args) < 1:
        return 1
    qual = process_fastq(args[0],PHRED_OFFSET)
    #print_average_quals(qual)

def print_average_quals(qualities):
    
    first = True
    for bin, avg_quality in sorted(qualities.items()):
        if first:
            print ",".join([str(i) for i in range(len(avg_quality)+1)])
            first = False
        avg_quality.insert(0,bin)
        print ",".join([str(i) for i in avg_quality])
        
def process_fastq(fastq_file, phred_offset=33):
    
    qual = {}
    count = {}
    out_handle = {}
    fh = miseq_data.FastQParser(fastq_file)
    outroot, ext = os.path.splitext(fastq_file)
    
    for record in fh:
        quality = record[3].strip()
        bin = int(round(quality_average(quality,phred_offset)))
        if bin not in out_handle:
            outfile = "%s.%d%s" % (outroot,bin,ext)
            out_handle[bin] = miseq_data.FastQWriter(outfile)
#            qual[bin] = [0. for i in range(len(quality))]
#            count[bin] = 0
#        for i,q in enumerate(quality):
#            qual[bin][i] += ord(q)
#        count[bin] += 1
        out_handle[bin].write(record)
        
#    for bin, sum_quality in qual.items():
#        for i,q in enumerate(sum_quality):
#            q -= phred_offset*count[bin]
#            qual[bin][i] = q/count[bin]
    
    for bin,oh in out_handle.items():
        oh.close()
          
    return qual


def quality_average(quality_str, phred_offset):
    sum = 0.
    for qual in quality_str:
        sum += ord(qual)
    sum -= len(quality_str)*phred_offset
    return sum/len(quality_str)
    
if __name__ == "__main__":
    sys.exit(main())
    