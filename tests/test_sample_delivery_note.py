import os
import unittest
import ConfigParser
from scilifelab.report.rl import make_example_note
from scilifelab.db.statusdb import SampleRunMetricsConnection

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

## Cutoffs
cutoffs = {
    "phix_err_cutoff" : 2.0,
    

    }


class TestSampleDeliveryNote(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(os.path.join(os.getenv("HOME"), "dbcon.ini")):
            self.url = None
            self.user = None
            self.pw = None
            self.examples = {}
        else:
            config = ConfigParser.ConfigParser()
            config.readfp(open(os.path.join(os.getenv("HOME"), "dbcon.ini")))
            self.url = config.get("couchdb", "url")
            self.user = config.get("couchdb", "username")
            self.pw = config.get("couchdb", "password")
            self.examples = {"sample":config.get("examples", "sample"),
                             "flowcell":config.get("examples", "flowcell"),
                             "project":config.get("examples", "project")}



    def test_1_make_example_note(self):
        make_example_note(os.path.join(filedir, "test.pdf"))


    def test_2_make_note(self):
        
