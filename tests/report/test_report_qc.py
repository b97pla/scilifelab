import os
import unittest
import logbook
from scilifelab.report.qc import fastq_screen, application_qc

from ..classes import has_couchdb_installation

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestReportQCFunctions(unittest.TestCase):
    def setUp(self):
        self.user = "user"
        self.pw = "pw"
        self.url = "localhost"
        self.examples = {"sample":"P001_101",
                         "flowcell":"120924_SN0002_0003_AC003CCCXX",
                         "project":"J.Doe_00_01"}

    def test_application_qc(self):
        """Test application specific qc"""
        data = application_qc(application="reseq", project_id=self.examples["project"], flowcell=self.examples["flowcell"].split("_")[-1], username=self.user, password=self.pw, projectdb="projects-test", sampledb="samples-test", url=self.url)
        tab = [x for x in data['stdout'].getvalue().split("\n") if x.startswith("P001")]
        self.assertEqual(len(tab), 2)

    def test_fqscreen(self):
        """Test fastq screen summary"""
        data = fastq_screen(project_id=self.examples["project"], flowcell=self.examples["flowcell"].split("_")[-1], username=self.user, password=self.pw, dbname="samples-test", url=self.url)
        self.assertEqual(len(data['stdout'].getvalue().split()), 2)
