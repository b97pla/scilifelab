"""
A module for handling the files and directories from an illumina run
"""

import os
import glob
from scilifelab.illumina.index_definitions import BASIC_LOOKUP
from scilifelab.utils.string import hamming_distance
 
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
        self.samplesheet = get_samplesheet(self.base)
        
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

    
    
        