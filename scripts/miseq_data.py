#!/usr/bin/env python
import os
import sys
import glob
import operator
import csv
import re
from optparse import OptionParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from bcbio.solexa.run_configuration import IlluminaConfiguration
from scilifelab.miseq import MiSeqRun

def main(run_dir):
    runobj = MiSeqRun(run_dir)
    runobj._split_fastq()
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--samplesheet", dest="samplesheet", default=None)
    options, args = parser.parse_args()
    main(args[0])
