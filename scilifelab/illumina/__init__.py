"""
A module for handling the files and directories from an illumina run
"""

import os
import glob
from scilifelab.illumina.index_definitions import BASIC_LOOKUP
from scilifelab.utils.string import hamming_distance
from scilifelab.bcbio.flowcell import Flowcell
from scilifelab.bcbio.qc import FlowcellRunMetricsParser
 
def map_index_name(index, mismatch=0):
    """Map the index sequences to the known names, if possible. Requires the samplesheet module.
    """
    
    names = []
    for name, sequence in BASIC_LOOKUP.items():
        sequence = sequence.replace('-','')
        try:
            if index != sequence and hamming_distance(index,sequence) > mismatch:
                continue
            names.append(name)
        except:
            pass
        
    return names

class IlluminaRun(object):
    
    def __init__(self, run_dir, samplesheet=None):
        self._run_dir = os.path.normpath(run_dir)
        assert os.path.exists(self._run_dir), "The path %s is invalid" % self._run_dir
        
        # Parse the run parameters
        parser = FlowcellRunMetricsParser(self._run_dir)
        self.run_config = parser.parseRunParameters()
        self.run_info = parser.parseRunInfo()
        
        self.samplesheet_file = samplesheet or IlluminaRun.get_samplesheet(self._run_dir)
    
    def get_base_mask(self):
        """Create the appropriate base mask based on the run setup
        """
        read_mask = []
        for read in self.run_info.get('Reads',[]):
            cycles = read.get('NumCycles',0)
            if read.get('IsIndexedRead','N') == 'Y':
                if cycles == "7":
                    cycles = "6N"
                read_mask.append("I{}".format(str(cycles)))
            else:
                read_mask.append("Y{}".format(str(cycles)))
                
        return ",".join(read_mask)
    
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
            elif ssheet is None:
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
        glob_pattern = os.path.join(self.get_sequence_dir_pattern(),"Project_*")
        for dir in glob.glob(glob_pattern):
            if project is not None:
                pname = os.path.basename(dir).split("_",1)[1].replace('__','.')
                if pname != project.replace('__','.'):
                    continue
            dirs.append(dir)    
        
        return dirs

    def get_sequence_dir_pattern(self):
        """Return the path pattern to the top directories containing the sequence reads
        """
        return os.path.join(self._run_dir,"Unaligned*")
    
    def get_sequence_dir(self):
        """Returns a list of the paths to the top directories containing the sequence reads
        """
        dirs = []
        for dir in glob.glob(self.get_sequence_dir_pattern()):
            if os.path.isdir(dir):
                dirs.append(dir)
        return dirs
        
    def get_unmatched_dir_pattern(self):
        """Returns the path pattern to the folder containing undetermined index reads
        """
        return os.path.join(self.get_sequence_dir_pattern(),"Undetermined_indices")
    
    def get_unmatched_dir(self):
        """Returns a list of the paths to the folders containing undetermined index reads
        """
        dirs = []
        for dir in glob.glob(self.get_unmatched_dir_pattern()):
            if os.path.isdir(dir):
                dirs.append(dir)
        return dirs
    
    def get_unmatched_reads(self, lanes=range(1,9)):
        """Return a list of fastq files with unmatched reads for each lane specified
        """
        
        reads = []
        for lane in lanes:
            fq_pattern = os.path.join(self.get_unmatched_dir_pattern(),"Sample_lane{:d}".format(lane),"lane{l:d}_Undetermined_L00{l:d}_R[12]_*.fastq.gz".format(l=lane))
            reads.append(glob.glob(fq_pattern))
        
        return reads
    
    def get_basecall_stats(self):
        """Return the path to the Basecall_stats_FCID directory
        """
        basecall_stats_dir_pattern = os.path.join(self.get_sequence_dir_pattern(),"Basecall_Stats_*")
        return glob.glob(basecall_stats_dir_pattern)
        
    def parse_directory(self):
        """Traverse a CASAVA 1.8+ generated directory structure and return a dictionary
        """ 
        raise NotImplementedError("This method is not yet implemented")
    
        projects = []
        
        # Create a Flowcell object
        fc = Flowcell()
        fc.filename = self.samplesheet_file
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
    
   
    
        
