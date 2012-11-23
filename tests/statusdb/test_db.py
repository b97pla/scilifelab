import os
import unittest
import ConfigParser
import logbook
from scilifelab.db.statusdb import SampleRunMetricsConnection

from ..classes import has_couchdb_installation

filedir = os.path.abspath(__file__)
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestDbConnection(unittest.TestCase):
    def setUp(self):
        self.user = "user"
        self.pw = "pw"
        self.url = "localhost"
        self.examples = {"sample":"1_120924_AC003CCCXX_TGACCA",
                         "flowcell":"AC003CCCXX",
                         "project":"J.Doe_00_01"}

    def test_connection(self):
        """Test database connection"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        self.assertEqual(sample_con.url_string, "http://{}:5984".format(self.url))

    def test_get_flowcell(self):
        """Test getting a flowcell for a given sample"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        fc = sample_con.get_entry(self.examples["sample"], "flowcell")
        self.assertEqual(str(fc), self.examples["flowcell"])

    def test_get_sample_ids(self):
        """Test getting sample ids given flowcell and sample_prj"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(sample_ids)))
        self.assertEqual(len(sample_ids), 5)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info( "Number of samples after subsetting: " + str(len(sample_ids)))
        self.assertEqual(len(sample_ids), 2)

    def test_get_samples(self):
        """Test getting samples given flowcell and sample_prj."""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)

        samples = sample_con.get_samples(fc_id=self.examples["flowcell"])
        LOG.info("Selecting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 5)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info("Selecting on flowcell, subsetting on project: " + str(len(samples)))
        self.assertEqual(len(samples), 2)

        samples = sample_con.get_samples(sample_prj=self.examples["project"])
        LOG.info("Selecting on project: " + str(len(samples)))
        self.assertEqual(len(samples), 3)
        samples = sample_con.get_samples(sample_prj=self.examples["project"], fc_id=self.examples["flowcell"])
        LOG.info("Selecting on project, subsetting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 2)

    def test_get_samples_wrong_info(self):
        """Test getting samples when either flowcell or project id information is wrong"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)

        samples = sample_con.get_samples(sample_prj="bogusproject", fc_id=self.examples["flowcell"])
        LOG.info("Selecting on bogus project, subsetting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 0)
        
                
    def test_get_project_sample_ids(self):
        """Test getting project sample ids"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(sample_prj=self.examples["project"])
        sample_names = [sample_con.db.get(x)["name"] for x in sample_ids]
        self.assertEqual(set(sample_names) , set(['1_120924_AC003CCCXX_TGACCA', '2_120924_AC003CCCXX_ACAGTG', '1_121015_BB002BBBXX_TGACCA']))
        
