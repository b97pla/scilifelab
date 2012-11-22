import os
import unittest
import logbook
import json
from scilifelab.db.statusdb import flowcell_run_metrics, project_summary, sample_run_metrics

filedir = os.path.abspath(__file__)
LOG = logbook.Logger(__name__)

class TestDocuments(unittest.TestCase):
    def test_doc(self):
        fc = flowcell_run_metrics(fc_name="Test", fc_date="121101")
        self.assertEqual(fc["name"], "121101_Test")
        self.assertEqual(fc["entity_type"], "flowcell_run_metrics")
        ps = project_summary()
        self.assertEqual(ps["entity_type"], "project_summary")
        kw = {'name':'test', 'barcode_name':'noeu'}
        srm = sample_run_metrics(**kw)
        self.assertEqual(srm["name"], "None_None_None_None")
