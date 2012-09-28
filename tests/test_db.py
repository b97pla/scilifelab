import os
import unittest
import ConfigParser
from scilifelab.db.statusdb import SampleRunMetricsConnection

filedir = os.path.abspath(__file__)

class TestDbConnection(unittest.TestCase):
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

    def test_1_connection(self):
        """Test database connection"""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        self.assertEqual(sample_con.url, "http://{}:5984".format(self.url))

    def test_2_get_flowcell(self):
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc = sample_con.get_entry(self.examples["sample"], "flowcell")
        self.assertEqual(str(fc), self.examples["flowcell"])

    def test_3_get_sample_ids(self):
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"])
        print len(sample_ids)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        print len(sample_ids)

    def test_4_get_samples(self):
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"])
        print len(samples)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        print len(samples)
                

