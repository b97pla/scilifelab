"""
A module for handling HiSeq-specificfiles and folders
"""
import os
import glob
import csv
import scilifelab.illumina as illumina

class HiSeqRun(illumina.IlluminaRun):
    
    def __init__(self, base, samplesheet=None):
        self.base = base
        if samplesheet is None:
            samplesheet = illumina.IlluminaRun.get_samplesheet(self.base)
        self.samplesheet = samplesheet

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
    
    def _unmatched_dir(self):
        """Returns the path to the folder containing undetermined index reads
        """
        return os.path.join(self.base,"Unaligned","Undetermined_indices")
    
    def get_unmatched_reads(self, lanes=range(1,9)):
        """Return a list of fastq files with unmatched reads for each lane specified
        """
        
        reads = []
        for lane in lanes:
            fq_pattern = os.path.join(self._unmatched_dir(),"Sample_lane{:d}".format(lane),"lane{l:d}_Undetermined_L00{l:d}_R[12]_*.fastq.gz".format(l=lane))
            reads.append(glob.glob(fq_pattern))
        
        return reads
    
             