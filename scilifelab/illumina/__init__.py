"""
A module for handling the files and directories from an illumina run
"""

import os
import glob
from scilifelab.illumina.index_definitions import BASIC_LOOKUP
from scilifelab.utils.string import hamming_distance
from scilifelab.bcbio.flowcell import Flowcell
 
def map_index_name(index, mismatch=0):
    """Map the index sequences to the known names, if possible. Requires the samplesheet module.
    """
    
    names = []
    for name, sequence in BASIC_LOOKUP.items():
        try:
            if index != sequence and hamming_distance(index,sequence) > mismatch:
                continue
            names.append(name)
        except:
            pass
        
    return names

class IlluminaRun():
    
    def __init__(self, base):
        self.base = base
        self.samplesheet = IlluminaRun.get_samplesheet(self.base)
        
    @staticmethod
    def get_samplesheet(base_dir):
        """Get the samplesheet from a flowcell directory, returning firstly [FCID].csv and secondly SampleSheet.csv
        """
        pattern = os.path.join(base_dir,"*.csv")
        ssheet = None
        for f in glob.glob(pattern):
            if not os.path.isfile(f):
                continue
            name, _ = os.path.splitext(os.path.basename(f))
            if base_dir.endswith(name):
                return f
            if name == "SampleSheet":
                ssheet = f
        return ssheet

    @staticmethod    
    def get_flowcell(root_dir, fc_dir=None):
        """Get the flowcells matching the (potentially partial) flowcell pattern
        """
        if fc_dir is None:
            fc_dir = ""
        
        # Match the pattern
        dir_pattern = os.path.join(root_dir,"*{}".format(fc_dir))
        dir_match = glob.glob(dir_pattern)
        
        # If the match is ambiguous but one exactly matches the pattern, return that
        for match in dir_match:
            if os.path.basename(match) == fc_dir:
                return [match]
        return dir_match

    def get_project_dir(self, project=None):
        """Return the path to the top directory containing the sample sequence reads for each project
        or, optionally, filtering by project
        """
        dirs = []
        glob_pattern = os.path.join(self.get_sequence_dir(),"Project_*")
        for dir in glob.glob(glob_pattern):
            if project is not None:
                pname = os.path.basename(dir).split("_",1)[1].replace('__','.')
                if pname != project.replace('__','.'):
                    continue
            dirs.append(dir)    
        
        return dirs

    def get_sequence_dir(self):
        """Return the path to the top directory containing the sequence reads
        """
        return os.path.join(self.base,"Unaligned")
        
    def get_unmatched_dir(self):
        """Returns the path to the folder containing undetermined index reads
        """
        return os.path.join(self.get_sequence_dir(),"Undetermined_indices")
    
    def get_unmatched_reads(self, lanes=range(1,9)):
        """Return a list of fastq files with unmatched reads for each lane specified
        """
        
        reads = []
        for lane in lanes:
            fq_pattern = os.path.join(self.get_unmatched_dir(),"Sample_lane{:d}".format(lane),"lane{l:d}_Undetermined_L00{l:d}_R[12]_*.fastq.gz".format(l=lane))
            reads.append(glob.glob(fq_pattern))
        
        return reads
    
    def get_basecall_stats(self):
        """Return the path to the Basecall_stats_FCID directory
        """
        basecall_stats_dir_pattern = os.path.join(self.get_sequence_dir(),"Basecall_Stats_*")
        basecall_stats_dir = glob.glob(basecall_stats_dir_pattern)
        if len(basecall_stats_dir) > 1:
            raise ValueError("ambiguous Basecall_Stats directories: {}".format("; ".join(basecall_stats_dir)))
        if len(basecall_stats_dir) == 1:
            return basecall_stats_dir[0]
        return None
        
    def parse_directory(self):
        """Traverse a CASAVA 1.8+ generated directory structure and return a dictionary
        """ 
        raise NotImplementedError("This method is not yet implemented")
    
        projects = []
        
        # Create a Flowcell object
        fc = Flowcell()
        fc.filename = self.samplesheet
        fc.fc_id()
        
        unaligned_dir = self.get_sequence_dir()
        basecall_stats_dir = self.get_basecall_stats()
        
        project_dir_pattern = os.path.join(unaligned_dir,"Project_*")
        for project_dir in glob.glob(project_dir_pattern):
            project_samples = []
            sample_dir_pattern = os.path.join(project_dir,"Sample_*")
            for sample_dir in glob.glob(sample_dir_pattern):
                fastq_file_pattern = os.path.join(sample_dir,"*.fastq.gz")
                samplesheet_pattern = os.path.join(sample_dir,"*.csv")
                fastq_files = [os.path.basename(file) for file in glob.glob(fastq_file_pattern)]
                samplesheet = glob.glob(samplesheet_pattern)
                assert len(samplesheet) == 1, "ERROR: Could not unambiguously locate samplesheet in %s" % sample_dir
                sample_name = sample_dir.replace(sample_dir_pattern[0:-1],'')
                project_samples.append({'sample_dir': os.path.relpath(sample_dir,project_dir), 'sample_name': sample_name, 'files': fastq_files, 'samplesheet': os.path.basename(samplesheet[0])})
            project_name = project_dir.replace(project_dir_pattern[0:-1],'')
            projects.append({'project_dir': os.path.relpath(project_dir,unaligned_dir), 'project_name': project_name, 'samples': project_samples})
        
        return {'fc_dir': fc_dir, 'fc_name': fc_name, 'fc_date': fc_date, 'data_dir': os.path.relpath(unaligned_dir,fc_dir), 'basecall_stats_dir': basecall_stats_dir, 'projects': projects}
    
   
    
        