import os
import unittest
import ConfigParser

import pandas as pd
from pandas.core.format import set_eng_float_format
from scilifelab.db.statusdb import SampleRunMetricsConnection
from test_default import SciLifeTest

filedir = os.path.abspath(__file__)

class TestQCData(SciLifeTest):
    def setUp(self):
        super(TestQCData, self).setUp()
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
            self.examples = {"sample":config.get("examples", "sample"),
                             "flowcell":config.get("examples", "flowcell"),
                             "project":config.get("examples", "project")}

    def tearDown(self):
        pass

    def test_1_qc(self):
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        qcdata = s_con.get_qc_data(self.examples["project"], self.examples["flowcell"])
        qcdf = pd.DataFrame(qcdata)
        print qcdf
        set_eng_float_format(accuracy=1, use_eng_prefix=True)
        qcdf.ix["TOTAL_READS"] = qcdf.ix["TOTAL_READS"] 
        qcdft = qcdf.T
        print qcdft.to_string()
