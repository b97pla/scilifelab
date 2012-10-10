"""Demultiplex a CASAVA 1.8+ FastQ file based on the information in the fastq header
"""

import sys
import argparse
from scilifelab.utils.fastq_utils import FastQParser, is_read_pair

def demultiplex_fastq(index, fastq1, fastq2):

    filter = {'index': index}
    
    fp1 = FastQParser(fastq1,filter)
    if fastq2 is not None:
        fp2 = FastQParser(fastq2,filter)
    for r1 in fp1:
        if fastq2 is not None:
            r2 = fp2.next()
            assert is_read_pair(r1,r2), "Mismatching headers for expected read pair" 
            sys.stderr.write("{}\n".format("\n".join(r2)))
             
        sys.stdout.write("{}\n".format("\n".join(r1))) 
    
def main():
    
    parser = argparse.ArgumentParser(description="Demultiplex a CASAVA 1.8+ FastQ file based on the information in the fastq header")

    parser.add_argument('fastq1', action='store', 
                        help="FastQ file to demultiplex")
    parser.add_argument('-f','--fastq2', action='store', default=None, 
                        help="Optional paired FastQ file to demultiplex")
    parser.add_argument('index', action='store', 
                        help="Index sequence to demultiplex on")
    
    args = parser.parse_args()
    demultiplex_fastq(args.index,args.fastq1,args.fastq2)
      
if __name__ == "__main__":
    main()
        