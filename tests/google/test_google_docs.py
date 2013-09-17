
import os
import tempfile
import shutil
import unittest
from scilifelab.google import get_credentials
from scilifelab.report.gdocs_report import _demultiplex_spreadsheet, _run_setup
import scilifelab.google.google_docs as gdocs

class TestGoogleDocs(unittest.TestCase):
    
    def setUp(self): 
        # Assume credentials can be found in ~/.pm/pm.conf
        cfg = os.path.expanduser(os.path.join("~",".pm","pm.conf"))
        assert os.path.exists(cfg), "{} configuration file does not exist, no gdocs credentials could be found"
        cfile = None
        with open(cfg) as fh:
            for line in fh:
                if not line.strip() == "[gdocs]":
                    continue
                for line in fh:
                    if line.strip()[0] == "[":
                        break
                    key, val = [s.strip() for s in line.split("=")]
                    if key == "credentials_file":
                        cfile = val
                        break
        assert cfile is not None, "Could not determine gdocs credentials file"
        self.credentials_file = os.path.expanduser(cfile)
        self.workdir = tempfile.mkdtemp("test_google_docs_")
        
    def tearDown(self):
        shutil.rmtree(self.workdir)
    
    def test_get_credentials(self):
        """Parse credentials file"""
        cfile = "credentials.tmp"
        exp_credentials = "credentials string"
        with open(cfile,"w") as fh:
            fh.write(exp_credentials)
        
        obs_credentials = get_credentials(cfile)
        self.assertEqual(exp_credentials,
                         obs_credentials,
                         "Expected and observed credentials are not identical")
        
    def test__demultiplex_spreadsheet(self):
        """Generate correct demultiplex spreadsheet name"""
        
        dates = ["100101",
                 "100331",
                 "110401",
                 "110630",
                 "120701",
                 "120930",
                 "131001",
                 "131231"]
        
        exp_names = ["Demultiplex Counts 2010 Q1",
                     "Demultiplex Counts 2010 Q1",
                     "Demultiplex Counts 2011 Q2",
                     "Demultiplex Counts 2011 Q2",
                     "Demultiplex Counts 2012 Q3",
                     "Demultiplex Counts 2012 Q3",
                     "Demultiplex Counts 2013 Q4",
                     "Demultiplex Counts 2013 Q4"]
        
        for i in range(len(dates)):
            obs_name = _demultiplex_spreadsheet(dates[i])
            self.assertEqual(exp_names[i],
                             obs_name,
                             "Generated name '{}' does not match expected name '{}'".format(obs_name,exp_names[i]))
    
    def test_Spreadsheet(self):
        """Instantiate and authenticate a SpreadSheet object"""
        credentials = get_credentials(self.credentials_file)
        obs_type = type(gdocs.SpreadSheet(credentials))
        exp_type = gdocs.SpreadSheet
        self.assertIs(obs_type,
                      exp_type,
                      "Did not instantiate an object of the correct type")

    def test_Document(self):
        """Instantiate and authenticate a Document object"""
        credentials = get_credentials(self.credentials_file)
        obs_type = type(gdocs.Document(credentials))
        exp_type = gdocs.Document
        self.assertIs(obs_type,
                      exp_type,
                      "Did not instantiate an object of the correct type")
         
    def test__run_setup(self):
        """Test getting run setup"""
        
        test_setups = {'1x50': [{'NumCycles': '50'}, {'IsIndexedRead': 'Y'}],
                       '2x50': [{'NumCycles': '50', 'IsIndexedRead': 'N'}, {'IsIndexedRead': 'Y'}, {'IsIndexedRead': 'Y'}, {'NumCycles': '50'}],
                       '100,50': [{'NumCycles': '100', 'IsIndexedRead': 'N'}, {'IsIndexedRead': 'Y'}, {'IsIndexedRead': 'Y'}, {'NumCycles': '50'}],
                       '2x?': [{'IsIndexedRead': 'N'}, {'IsIndexedRead': 'Y'}, {'IsIndexedRead': 'Y'}, {}]}
        
        for exp_setup, reads in test_setups.items():
            obs_setup = _run_setup(reads)
            self.assertEqual(exp_setup,
                             obs_setup,
                             "Observed run setup '{}' does not match expected setup '{}'".format(obs_setup,exp_setup))
