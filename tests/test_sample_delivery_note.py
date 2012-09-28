import os
import unittest
import ConfigParser
from scilifelab.report.rl import make_example_note
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

## Cutoffs
cutoffs = {
    "phix_err_cutoff" : 2.0,
    }

## parameters
parameters = {
    "customer_reference": None,
    "FC_id" : None,
    "scilifelab_name" : None,
    "rounded_read_count" : None,
    "phix_error_rate" : None,
    "avg_quality_score" : None,
    "ordered_amount" : None,
    "success" : None,
    "project_name" : None,
}

## key mapping from sample_run_metrics to parameter keys
srm_to_parameter = {"project_name":"sample_prj", "FC_id":"flowcell", 
                    "scilifelab_name":"barcode_name"}
## mapping project_summary to parameter keys
ps_to_parameter = {"":""}

## Data dicts
error_rates = {}
qvs = {}



class TestSampleDeliveryNote(unittest.TestCase):
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
            self.examples = {"sample":config.get("examples", "sample"),
                             "flowcell":config.get("examples", "flowcell"),
                             "project":config.get("examples", "project")}



    def test_1_make_example_note(self):
        """Make example note"""
        make_example_note(os.path.join(filedir, "test.pdf"))


    def test_2_make_note(self):
        """Make a note subset by example flowcell and project"""
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc_con = FlowcellRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = s_con.get_samples(self.examples["flowcell"], self.examples["project"])
        fc = "{}_{}".format(samples[0]["date"], samples[0]["flowcell"])
        fc_id = fc_con.get_phix_error_rate(str(fc))
        # for s in samples:
        #     s_param = parameters
        #     print s_param
        #     s_param.update({s})
        #     print s_param
