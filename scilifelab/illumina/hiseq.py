"""
A module for handling HiSeq-specificfiles and folders
"""
import csv

class HiSeq(IlluminaRun):
    pass

    def parse_samplesheet(self):
        """Parse a .csv samplesheet and return a list of dictionaries with
        elements corresponding to rows of the samplesheet and keys corresponding
        to the columns in the header
        """
        entries = []
        with open(self.samplesheet) as fh:
            csvread = csv.DictReader(fh, dialect='excel')
            entries = [row for row in csvread]
        
        return entries
    