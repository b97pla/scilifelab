import os
import unittest
import logbook
from scilifelab.report import sequencing_success
from scilifelab.report.rl import make_example_sample_note

from ..classes import has_couchdb_installation

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestSampleDeliveryNote(unittest.TestCase):
    def setUp(self):
        self.user = "user"
        self.pw = "pw"
        self.examples = {"sample":"P001_101",
                         "flowcell":"120924_SN0002_0003_AC003CCCXX",
                         "project":"J.Doe_00_01"}

    def test_make_example_note(self):
        """Make example note"""
        success = True
        try:
            make_example_sample_note(os.path.join(filedir, "test.pdf"))
        except:
            success = False
        self.assertTrue(success)

class TestSampleDeliveryFunctions(unittest.TestCase):
    def test_sequencing_success(self):
        """Make sure sequencing success returns correct assessment"""
        cutoffs = {'phix_err_cutoff':2.0}
        # Successful run
        msg = sequencing_success({'phix_error_rate':1.3, 'rounded_read_count':50.3, 'ordered_amount':35.2}, cutoffs)
        self.assertEqual(msg, "Successful run.")
        # Failed run on account of # sequences
        msg = sequencing_success({'phix_error_rate':1.3, 'rounded_read_count':11.3, 'ordered_amount':35.2}, cutoffs)
        self.assertEqual(msg, "The yield may be lower than expected.")
        # Failed run on account of phix error rate
        msg = sequencing_success({'phix_error_rate':2.3, 'rounded_read_count':50.3, 'ordered_amount':35.2}, cutoffs)
        self.assertEqual(msg, "High average error rate.")
        # Failed run on account of both phix error rate and # sequences
        msg = sequencing_success({'phix_error_rate':2.3, 'rounded_read_count':20.3, 'ordered_amount':35.2}, cutoffs)
        self.assertEqual(msg, "High average error rate.The yield may be lower than expected.")


        
