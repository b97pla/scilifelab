
import unittest
import shutil
import tempfile
import os
import random
import string
import bcbio.utils as utils
import tests.generate_test_data as td
from scilifelab.illumina.hiseq import HiSeqRun

class TestHiSeqRun(unittest.TestCase):
    
    def setUp(self): 
        self.rootdir = tempfile.mkdtemp(prefix="test_illumina_hiseq_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
    
    def test_parse_samplesheet(self):
        """Write and parse a csv-file
        """
        
        # Assert non-existing file raises exception
        with self.assertRaises(IOError):
            HiSeqRun.parse_samplesheet(os.path.join(self.rootdir,'non-existing-samplesheet'))
            
        # Write a csv file with some bogus values
        sdata = td.generate_samplesheet_data()
        samplesheet = os.path.join(self.rootdir,'SampleSheet.csv')
        HiSeqRun.write_samplesheet(sdata,samplesheet)
        
        # Assert that the written data corresponds to the generated data
        with open(samplesheet) as fh:
            # Assert that header is correct
            self.assertListEqual(HiSeqRun._samplesheet_header(),
                                 fh.next().strip().split(","),
                                 "Written header does not match expected header")
            for entry in sdata:
                # Assert that all rows have the correct values in the correct columns
                self.assertListEqual([str(e) for e in entry],
                                     fh.next().strip().split(","),
                                     "Written data row does not match entry in generated samplesheet")
            
            # Assert that all rows from samplesheet has been consumed
            with self.assertRaises(StopIteration):
                fh.next()
        
        # Assert that the parsed data matches the generated data
        data = HiSeqRun.parse_samplesheet(samplesheet)
        self.assertEqual(len(sdata),
                         len(data),
                         "Number of parsed entries does not match number of generated entries")
        for d in data:
            self.assertListEqual([str(e) for e in sdata.pop(0)],
                                 [d[col] for col in HiSeqRun._samplesheet_header()],
                                 "Parsed data row does not match entry in generated samplesheet")
        
        