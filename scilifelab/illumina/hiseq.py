"""
A module for handling HiSeq-specificfiles and folders
"""
import csv
import scilifelab.illumina as illumina

class HiSeqRun(illumina.IlluminaRun):
    pass

    @staticmethod
    def _samplesheet_header():
        """Return a list of columns in the HiSeq samplesheet
        """
        return ["FCID",
                "Lane",
                "SampleID",
                "SampleRef",
                "Index",
                "Description",
                "Control",
                "Recipe",
                "Operator",
                "SampleProject"] 

    @staticmethod
    def parse_samplesheet(samplesheet):
        """Parse a .csv samplesheet and return a list of dictionaries with
        elements corresponding to rows of the samplesheet and keys corresponding
        to the columns in the header
        """
        entries = []
        with open(samplesheet) as fh:
            csvr = csv.DictReader(fh, dialect='excel')
            entries = [row for row in csvr]
        
        return entries
    
    @staticmethod
    def write_samplesheet(sdata, samplesheet):
        """Write a .csv samplesheet from a list of entries
        """
        with open(samplesheet,"w") as outh:
            csvw = csv.writer(outh)
            csvw.writerow(HiSeqRun._samplesheet_header())
            csvw.writerows(sdata)
        return samplesheet
    