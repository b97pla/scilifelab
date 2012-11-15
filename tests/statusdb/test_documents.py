import os
import unittest
import logbook
import json
from scilifelab.db.statusdb import FlowcellRunMetrics

filedir = os.path.abspath(__file__)
LOG = logbook.Logger(__name__)

class TestDocuments(unittest.TestCase):
    def test_doc(self):
        fc = FlowcellRunMetrics(fc_name="Test", fc_date="121101")
        print fc
        print dir(fc)
        print json.dumps(fc)
        print fc["entity_type"]
        print repr(fc)
