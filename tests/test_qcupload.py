import os
import unittest
import ConfigParser
from scilifelab.db.statusdb import SampleRunMetricsConnection
from scilifelab.bcbio.qc import FlowcellRunMetrics, SampleRunMetrics

filedir = os.path.abspath(__file__)

class TestQCUpload(unittest.TestCase):
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
            self.demuxstats = config.get("statusdb", "demuxstats")
            self.sample_kw = {"path":config.get("statusdb", "fcdir"),
                             "flowcell":config.get("statusdb", "fc_name"),
                             "date":config.get("statusdb", "date"),
                             "lane":config.get("statusdb", "lane"),
                             "barcode_name":config.get("statusdb", "name"),
                             "sample_prj":config.get("statusdb", "project"),
                             "barcode_id":config.get("statusdb", "barcode_id"),
                             "sequence":config.get("statusdb", "sequence")}
            self.fc_kw = {"path":config.get("statusdb", "fcdir"),
                          "fc_date":config.get("statusdb", "date"),
                          "fc_name":config.get("statusdb", "fc_name")}
            self.examples = {"sample":config.get("examples", "sample"),
                             "flowcell":config.get("examples", "flowcell"),
                             "project":config.get("examples", "project")}

    def test_1_demuxstats(self):
        obj = FlowcellRunMetrics(**self.fc_kw)
        metrics = obj.parse_demultiplex_stats_htm()
        print metrics["Barcode_lane_statistics"][0]
            
