"""Script for extracting and analyzing indexes found in the undetermined index read fraction
"""

import argparse
import os
import glob

from bcbio.pipeline.config_loader import load_config
from scilifelab.illumina import IlluminaRun

def analyze_indexes(config_file, flowcell, samplesheet=None, analysis_dir=None, brief=False):
    """High-level method for analyzing undetermined indexes
    """
    # Parse the config
    config = {}
    config = load_config(config_file)
    
    # Get the analysis dir from the config unless it was supplied
    if analysis_dir is None:
        analysis_dir = config.get("analysis",{}).get("base_dir",None)
        assert analysis_dir is not None, "The analysis directory could not be read from configuration"
    
    # Assert that the flowcell directory exists
    fc_dir = IlluminaRun._get_flowcell(analysis_dir,flowcell)
    assert len(fc_dir) > 0, "Could not locate flowcell *{} in analysis directory {}".format(flowcell,analysis_dir)
    assert len(fc_dir) == 1, "Ambiguous matches for flowcell *{} in analysis directory {}".format(flowcell,analysis_dir)
    fc_dir = fc_dir[0]

    # If not specified, locate the samplesheet in the analysis directory
    if samplesheet is None:
        samplesheet = IlluminaRun._get_samplesheet(fc_dir)
        
    

    
def main():
    """Argument parsing
    """
    
    parser = argparse.ArgumentParser(description="""Analyze the top N indexes among the reads with undetermined indexes and 
    flag any obvious problems detected""")

    parser.add_argument('-b','--brief', action='store_true', default=False, 
                        help="display only a brief summary of the status status")
    parser.add_argument('-a','--analysis-dir', dest='analysis_dir', action='store', default=None, 
                        help="path to the folder containing project analysis data. If specified, overrides the configuration value")
    parser.add_argument('-s','--samplesheet', action='store', default=None,
                        help=".csv samplesheet for the flowcell to be analyzed. If specified, overrides the samplesheet in default location")
    parser.add_argument('config', action='store', 
                        help="configuration file with paths, credentials etc.")
    parser.add_argument('flowcell', action='store', 
                        help="flowcell to analyze undetermined indexes for")
    args = parser.parse_args()
    analyze_indexes(args.config, args.flowcell, args.samplesheet, args.analysis_dir, args.brief)
      
if __name__ == "__main__":
    main()
 
