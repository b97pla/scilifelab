"""
Test runinfo functions
"""
import os
import yaml
import glob
import re
from cement.core import handler
from test_default import PmTest
from pmtools.lib.flowcell import *

flowcell = "120829_SN0001_0001_AA001AAAXX"
fc_dir = os.path.join(os.path.curdir, "data", "analysis", flowcell)
runinfo = os.path.join(os.path.curdir, "data", "archive", flowcell, "run_info.yaml")

class PmFlowcellTest(PmTest):
    """Test flowcell object functionality"""
    def test_1_get_flowcell(self):
        """Read contents of runinfo file and generate flowcell object"""
        fc = Flowcell(runinfo)
        newfc=fc.subset("sample_prj", "J.Doe_00_01")
        self.eq(len(fc), 11)
        self.eq(len(newfc), 7)
        self.eq(fc.projects(), ['J.Doe_00_01', 'J.Doe_00_02'])
        self.eq(newfc.projects(), ['J.Doe_00_01'])
        self.eq(fc.path, os.path.abspath(os.path.dirname(runinfo)))
        self.eq(fc.filename, os.path.abspath(runinfo))

        
    def test_2_load_flowcell(self):
        """Create object and load run information"""
        self.app = self.make_app(argv = [])
        self.app.setup()
        fc = Flowcell()
        fc.load([os.path.join(self.app.config.get("analysis", "root"), flowcell),
                 os.path.join(self.app.config.get("archive", "root"), flowcell)])
        self.eq(len(fc), 11)
        self.eq(fc.path, os.path.join(self.app.config.get("archive", "root"), flowcell))
        self.eq(fc.filename, os.path.join(self.app.config.get("archive", "root"), flowcell, "run_info.yaml"))

    def test_3_getters(self):
        """Test getters"""
        fc = Flowcell(runinfo)
        self.eq(fc.barcodes('2'), [5,7,17,19])
        self.eq(fc.lanes(), ['1','2'])

    def test_4_barcode_mapping(self):
        """Test barcode id to name mappings"""
        fc = Flowcell(runinfo)
        self.eq(fc.barcode_id_to_name('1'), {1: 'P1_101F_index1', 2: 'P1_102F_index2', 3: 'P1_103_index3', 4: 'P1_104F_index4', 8: 'P1_105F_index5', 10: 'P1_106F_index6', 12: 'P1_107_index7'})
        self.eq(fc.barcode_name_to_id('1'), {'P1_102F_index2': 2, 'P1_107_index7': 12, 'P1_106F_index6': 10, 'P1_105F_index5': 8, 'P1_101F_index1': 1, 'P1_103_index3': 3, 'P1_104F_index4': 4})


    def test_5_glob_str(self):
        """Test construction of glob prefixes"""
        fc = Flowcell(runinfo)
        glob_pfx_str = fc.glob_pfx_str()
        self.app = self.make_app(argv = [])
        self.app.setup()
        glob_str = os.path.join(self.app.config.get("analysis", "root"), flowcell, glob_pfx_str[0])
        print glob_str
        res = re.search(glob_str)
        print res
                
    def test_6_get_files(self):
        """Test getting files"""
        fc = Flowcell(runinfo)
        #fc.get_files(fc_dir)
        fc.get_files(fc_dir, project="J.Doe_00_01")
