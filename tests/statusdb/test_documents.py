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
        print fc
        print dir(fc)
        print json.dumps(fc)
        print fc["entity_type"]
        print repr(fc)
        ps = project_summary()
        print ps
        print json.dumps(ps)
        kw = {'name':'test', 'barcode_name':'noeu'}
        srm = sample_run_metrics(**kw)

        print json.dumps(srm)
