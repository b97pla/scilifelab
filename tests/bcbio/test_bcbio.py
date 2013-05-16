import os
import tempfile
import shutil
import unittest
from ..data import data_files
from scilifelab.bcbio.qc import RunInfoParser

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

RunInfo = data_files["RunInfo.xml"]

class TestBcbioQC(unittest.TestCase):
    """Test for bcbio qc module"""
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_bcbio_qc_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
    
    def test_parse_runinfo(self):
        temp = tempfile.TemporaryFile(mode="w+t")
        temp.write(RunInfo)
        temp.seek(0)
        rip = RunInfoParser()
        res = rip.parse(temp)
        self.assertEqual(res["Id"], "120924_SN0002_0003_CC003CCCXX")
        self.assertEqual(res["Flowcell"], "CC003CCCXX")
        self.assertEqual(res["Instrument"], "SN0002")
        self.assertEqual(res["Date"], "120924")

