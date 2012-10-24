import os
import unittest
import ConfigParser
import logbook
from scilifelab.db.statusdb import SampleRunMetricsConnection

filedir = os.path.abspath(__file__)
LOG = logbook.Logger(__name__)

@unittest.skipIf(not os.path.exists(os.path.join(os.getenv("HOME"), "dbcon.ini")), "database tests require a database connection config file ~/dbcon.ini")
class TestDbConnection(unittest.TestCase):
    def setUp(self):
        config = ConfigParser.ConfigParser()
        config.readfp(open(os.path.join(os.getenv("HOME"), "dbcon.ini")))
        self.url = config.get("couchdb", "url")
        self.user = config.get("couchdb", "username")
        self.pw = config.get("couchdb", "password")
        self.examples = {"sample":config.get("examples", "sample"),
                         "flowcell":config.get("examples", "flowcell"),
                         "project":config.get("examples", "project")}

    def test_connection(self):
        """Test database connection"""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        self.assertEqual(sample_con.url_string, "http://{}:5984".format(self.url))

    def test_get_flowcell(self):
        """Test getting a flowcell for a given sample"""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc = sample_con.get_entry(self.examples["sample"], "flowcell")
        self.assertEqual(str(fc), self.examples["flowcell"])

    def test_get_sample_ids(self):
        """Test getting sample ids given flowcell and sample_prj"""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(sample_ids)))
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info( "Number of samples after subsetting: " + str(len(sample_ids)))

    def test_get_samples(self):
        """Test getting samples given flowcell and sample_prj."""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(samples)))
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info("Number of samples after subsetting: " + str(len(samples)))
                
    def test_get_project_sample_ids(self):
        """Test getting project sample ids"""
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(sample_prj=self.examples["project"])
        fc_sample_ids = ["sample1"]
        prj_sample_ids = ["sample1", "sample2"]
        ids = list(set(fc_sample_ids) | set(prj_sample_ids))
        self.assertEqual(ids, prj_sample_ids)
