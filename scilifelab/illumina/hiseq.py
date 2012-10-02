"""
A module for handling HiSeq-specificfiles and folders
"""
import csv
import scilifelab.illumina as illumina

class HiSeqRun(illumina.IlluminaRun):
    
    def __init__(self, base, samplesheet=None):
        self.base = base
        if samplesheet is None:
            samplesheet = illumina.IlluminaRun.get_samplesheet(self.base)
        self.samplesheet = samplesheet
        self.project_names = get_project_names(self.samplesheet)

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
    
    @staticmethod
    def get_project_names(samplesheet):
        """List the projects available in the samplesheet. Optionally filter by project name.
        """ 
        return sorted(list(set([e['SampleProject'].replace("__",".") for e in HiSeqRun.parse_samplesheet(samplesheet)])))
    
    @staticmethod
    def get_project_sample_ids(samplesheet, project):
        """Return the samples listed in the samplesheet for a project
        """
        ids = []
        for e in HiSeqRun.parse_samplesheet(samplesheet):
            if e['SampleProject'].replace('__','.') == project:
                ids.append(e['SampleID'])
        return ids
    