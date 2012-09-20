import glob
import os
import collections
import argparse

def main():
    
    parser = argparse.ArgumentParser(description="Query the status of flowcells, projects, samples that are organized "\
                                     "according to the CASAVA file structure")

    parser.add_argument('-f','--flowcell', action='store', default="*CXX", 
                        help="only display the status for the specified flowcell")
    parser.add_argument('-p','--project', action='store', default=None, 
                        help="display only the status for the specified project")
    parser.add_argument('-r','--archive-dir', dest='archive_dir', action='store', default="/proj/a2010002/archive", 
                        help="path to the folder containing flowcell data")
    parser.add_argument('-a','--analysis-dir', dest='analysis_dir', action='store', default="/proj/a2010002/nobackup/illumina", 
                        help="path to the folder containing project analysis data")
    parser.add_argument('run_folder', action='store', default=None, 
                        help="the full path to the run folder containing analysis output", nargs='+')
    
    args = parser.parse_args()
    status_query(args.archive_dir,args.analysis_dir,args.flowcell,args.project)
      
if __name__ == "__main__":
    main()

### ---- Tests ----

import unittest
import tempfile
import shutil
from mock import Mock

class TestMetadata(unittest.TestCase):
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="status_query_metadata_test")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
    