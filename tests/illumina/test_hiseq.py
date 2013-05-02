
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
        self.hiseq = HiSeqRun(self.rootdir)
        
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
            
        # Assert that filtering on lane returns expected output
        lanes = list(set([d["Lane"] for d in data]))
        obs_lane_data = HiSeqRun.parse_samplesheet(samplesheet,lane=lanes[-1])
        exp_lane_data = [d for d in data if str(d["Lane"]) == str(lanes[-1])]
        self.assertListEqual(sorted(obs_lane_data),
                             sorted(exp_lane_data),
                             "Parsed data row does not match entry in generated samplesheet")
            
    
    def test_get_project_names(self):
        """Get the projects from a samplesheet
        """
        # Assert that an empty file returns an empty list
        fh, ssheet = tempfile.mkstemp(dir=self.rootdir, suffix=".csv")
        os.close(fh)
        self.assertListEqual([],HiSeqRun.get_project_names(ssheet),
                             "The list of projects for an empty file is not empty")
        
        # Generate artificial samplesheet data
        data = td.generate_samplesheet_data()
        projects = {}
        for d in data:
            projects[d[-1]] = 1
        
        # Write the data to a samplesheet
        td._write_samplesheet(data,ssheet)
         
        # Assert that the list of projects returned is the same that we generated
        self.assertListEqual(sorted(projects.keys()),sorted(HiSeqRun.get_project_names(ssheet)),
                             "The list of projects does not match the original list")
        

    def test_get_project_sample_ids(self):
        """Test that getting the project samples from a samplesheet behaves as expected
        """
        
        # Generate artificial samplesheet data
        data = td.generate_samplesheet_data()
        fh, ssheet = tempfile.mkstemp(dir=self.rootdir, suffix=".csv")
        os.close(fh)
        td._write_samplesheet(data,ssheet)
         
        # Assert that getting samples for a non-existing project returns an empty list
        self.assertListEqual([],HiSeqRun.get_project_sample_ids(ssheet,td.generate_project()),
                             "Getting samples for a non-existing project returned unexpected output")
        
        # Iterate over the projects and assert that the returned samples are correct
        samples = {}
        for row in data:
            if row[9] not in samples:
                samples[row[9]] = []
            samples[row[9]].append(row[2])
        
        for proj, sample in samples.items():
            self.assertListEqual(sorted(sample),sorted(HiSeqRun.get_project_sample_ids(ssheet,proj)),
                                 "The returned list of samples did not match the original")
   
    