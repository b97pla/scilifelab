import os
import unittest
import logbook
from scilifelab.report.rl import make_example_project_note

from ..classes import has_couchdb_installation

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)


@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestProjectStatusNote(unittest.TestCase):
    def setUp(self):
        self.user = "user"
        self.pw = "pw"
        self.examples = {"sample":"P001_101",
                         "flowcell":"120924_SN0002_0003_AC003CCCXX",
                         "project":"J.Doe_00_01"}

    def test_make_example_project_note(self):
        """Make example project note"""
        success = True
        try:
            make_example_project_note(os.path.join(filedir, "project_test.pdf"))
        except:
            success = False
        self.assertTrue(success)

