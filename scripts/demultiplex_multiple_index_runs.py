"""Demultiplex a CASAVA 1.8+ FastQ file based on the information in the fastq header
"""

import sys
import argparse
from scilifelab.utils.fastq_utils import FastQParser, is_read_pair

def demultiplex_fastq(index, fastq1, fastq2):
    print "test\n";
    
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
        