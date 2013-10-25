""" A module for handling HiSeq-specific files and folders
"""
import os
import glob
import csv
import scilifelab.illumina as illumina


class HiSeqRun(illumina.IlluminaRun):
    def __init__(self, run_dir, samplesheet=None):
        illumina.IlluminaRun.__init__(self, run_dir, samplesheet)
        if self.samplesheet_file is not None:
            self.samplesheet = HiSeqSampleSheet(self.samplesheet_file)

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
    def parse_samplesheet(samplesheet, lane=None, sample_project=None, index=None):
        """Parse a .csv samplesheet and return a list of dictionaries with
        elements corresponding to rows of the samplesheet and keys
        corresponding to the columns in the header. Optionally filter by lane 
        and/or sample_project and/or index.
        """
        entries = []
        with open(samplesheet,"rU") as fh:
            csvr = csv.DictReader(fh, dialect='excel')
            entries = [row for row in csvr \
                       if (lane is None or row["Lane"] == lane) \
                       and (sample_project is None or row["SampleProject"] == sample_project) \
                       and (index is None or row["Index"] == index)]
        
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
        return sorted(list(set([e['SampleProject'].replace("__", ".") for e in HiSeqRun.parse_samplesheet(samplesheet)])))
    
    @staticmethod
    def get_project_sample_ids(samplesheet, project):
        """Return the samples listed in the samplesheet for a project
        """
        ids = []
        for e in HiSeqRun.parse_samplesheet(samplesheet):
            if e['SampleProject'].replace('__','.') == project:
                ids.append(e['SampleID'])
        return ids


class HiSeqSampleSheet(list):
    def __init__(self, samplesheet, lane=None, sample_project=None, index=None):
        self.header = ["FCID",
                       "Lane",
                       "SampleID",
                       "SampleRef",
                       "Index",
                       "Description",
                       "Control",
                       "Recipe",
                       "Operator",
                       "SampleProject"] 

        if isinstance(samplesheet, list):
            self.extend(samplesheet)

        else:
            self.samplesheet = samplesheet
            self._parse_sample_sheet(lane=None, sample_project=None, index=None)


    def _parse_sample_sheet(self, lane=None, sample_project=None, index=None):
        """Parse a .csv samplesheet and return a list of dictionaries with
        elements corresponding to rows of the samplesheet and keys
        corresponding to the columns in the header. Optionally filter by lane 
        and/or sample_project and/or index.
        """
        with open(self.samplesheet,"rU") as fh:
            csvr = csv.DictReader(fh, dialect='excel')
            for row in csvr:
                if (lane is None or row["Lane"] == lane) \
                and (sample_project is None or row["SampleProject"] == sample_project) \
                and (index is None or row["Index"] == index):
                    self.append(row)

    def write(self, samplesheet):
        """Write samplesheet to .csv file
        """
        with open(samplesheet, "w") as outh:
            csvw = csv.writer(outh)
            if len(self) > 0:
                csvw.writerow(self[0].keys())
            else:
                csvw.writerow(self.header)
            csvw.writerows([row.values() for row in self])
