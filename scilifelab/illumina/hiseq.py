"""
A module for handling HiSeq-specificfiles and folders
"""
import csv
from scilifelab.illumina import IlluminaRun

class HiSeqRun(IlluminaRun):
    pass

    @staticmethod
    def parse_samplesheet(samplesheet):
        """Parse a .csv samplesheet and return a list of dictionaries with
        elements corresponding to rows of the samplesheet and keys corresponding
        to the columns in the header
        """
        entries = []
        with open(samplesheet) as fh:
            csvread = csv.DictReader(fh, dialect='excel')
            entries = [row for row in csvread]
        
        return entries
    