"""
Test runinfo functions
"""
import os
import yaml
from cement.core import handler
from test_default import PmTest
from pmtools.lib.runinfo import *
from pmtools.lib.flowcell import *

runinfo = os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml")

class PmRuninfoTest(PmTest):
    def test_1_get_rows(self):
        """Given a runinfo object, retrieve rows based on different filtering criteria"""
        info = get_runinfo(runinfo, tab=True)
        subinfo = subset_runinfo(info, "sample_prj", "J.Doe_00_01")
        print dump_runinfo(subinfo)
        print dump_runinfo(info)

    def test_2_dump_runinfo(self):
        """Given a runinfo object, dump to cStringIO"""
        info = get_runinfo(runinfo, tab=True)
        print dump_runinfo(info)
        info = get_runinfo(runinfo, tab=False)
        print dump_runinfo(info, tab=False)


    def test_3_make_files(self):
        """From a runinfo table, make file names"""
        info = get_runinfo(runinfo, tab=True)
        res = find_files(info, "./data/analysis/120829_SN0001_0001_AA001AAAXX")
        print res

class PmFlowcellTest(PmTest):
    """Test flowcell object functionality"""
    def test_1_get_flowcell(self):
        """Read contents of runinfo file and generate flowcell object"""
        fc = Flowcell(runinfo)
        print fc
        print fc.as_yaml()
        newfc=fc.subset("sample_prj", "J.Doe_00_01")
        print newfc
