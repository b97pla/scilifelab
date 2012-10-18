import os
import unittest
import ConfigParser
import logbook
from scilifelab.db.statusdb import SampleRunMetricsConnection

filedir = os.path.abspath(__file__)
LOG = logbook.Logger(__name__)

class TestDbConnection(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(os.path.join(os.getenv("HOME"), "dbcon.ini")):
            self.url = None
            self.user = None
            self.pw = None
            self.examples = {}
            LOG.warning("No such file {}; will not run database connection tests".format(os.path.join(os.getenv("HOME"), "dbcon.ini")))
        else:
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
        if not self.examples:
            LOG.info("Not running test")
            return
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        self.assertEqual(sample_con.url_string, "http://{}:5984".format(self.url))

    def test_get_flowcell(self):
        """Test getting a flowcell for a given sample"""
        if not self.examples:
            LOG.info("Not running test")
            return
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc = sample_con.get_entry(self.examples["sample"], "flowcell")
        self.assertEqual(str(fc), self.examples["flowcell"])

    def test_get_sample_ids(self):
        """Test getting sample ids given flowcell and sample_prj"""
        if not self.examples:
            LOG.info("Not running test")
            return
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(sample_ids)))
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info( "Number of samples after subsetting: " + str(len(sample_ids)))

    def test_get_samples(self):
        """Test getting samples given flowcell and sample_prj."""
        if not self.examples:
            LOG.info("Not running test")
            return
        sample_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(samples)))
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info("Number of samples after subsetting: " + str(len(samples)))
                
