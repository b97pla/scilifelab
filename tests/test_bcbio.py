import os
import unittest
from pm.data import setup_data_files
from scilifelab.bcbio.qc import RunInfoParser

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

class TestBcbioQC(unittest.TestCase):
    """Test for bcbio qc module"""
    def setUp(self):
        setup_data_files()
    
    def test_1_parse_runinfo(self):
        infile = os.path.join("pm/data/archive/120924_SN0002_0003_CC003CCCXX/RunInfo.xml")
        rip = RunInfoParser()
        with open(infile) as fh:
            res = rip.parse(fh)
        self.assertEqual(res["Id"], "120924_SN0002_0003_CC003CCCXX")
        self.assertEqual(res["Flowcell"], "CC003CCCXX")
        self.assertEqual(res["Instrument"], "SN0002")
        self.assertEqual(res["Date"], "120924")

