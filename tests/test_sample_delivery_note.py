import os
import unittest
import ConfigParser
from scilifelab.report import sequencing_success, set_status
from scilifelab.report.rl import make_example_sample_note, make_note, sample_note_paragraphs, sample_note_headers
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectSummaryConnection

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

## Cutoffs
cutoffs = {
    "phix_err_cutoff" : 2.0,
    }

## parameters
parameters = {
    "project_name" : None,
    "customer_reference": None,
    "uppnex_project_id" : None,
    "ordered_amount" : None,
    "start_date" : None,
    "FC_id" : None,
    "scilifelab_name" : None,
    "customer_name" : None,
    "rounded_read_count" : None,
    "phix_error_rate" : None,
    "avg_quality_score" : None,
    "success" : None,
}

## key mapping from sample_run_metrics to parameter keys
srm_to_parameter = {"project_name":"sample_prj", "FC_id":"flowcell", 
                    "scilifelab_name":"barcode_name", "start_date":"date", "rounded_read_count":"bc_count"}
## mapping project_summary to parameter keys
ps_to_parameter = {"customer_reference":"customer_reference", "uppnex_project_id":"uppnex_id"}

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
        make_example_sample_note(os.path.join(filedir, "test.pdf"))

    def test_2_make_note(self):
        """Make a note subset by example flowcell and project"""
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc_con = FlowcellRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        paragraphs = sample_note_paragraphs()
        headers = sample_note_headers()
        samples = s_con.get_samples(self.examples["flowcell"], self.examples["project"])
        project = p_con.get_entry(self.examples["project"])
        samples = p_con.map_srm_to_name(self.examples["project"], fc_id=self.examples["flowcell"], use_bc_map=True)
        for k,v  in samples.items():
            s_param = parameters
            s = s_con.get_entry(k)
            s_param.update({key:s[srm_to_parameter[key]] for key in srm_to_parameter.keys()})
            fc = "{}_{}".format(s["date"], s["flowcell"])
            s_param["phix_error_rate"] = fc_con.get_phix_error_rate(str(fc), s["lane"])
            s_param['avg_quality_score'] = s_con.calc_avg_qv(s["name"])
            s_param['rounded_read_count'] = round(float(s_param['rounded_read_count'])/1e6,1) if s_param['rounded_read_count'] else None
            s_param['customer_name'] = project['samples'][v["sample"]].get('customer_name', None)

            if project:
                s_param['ordered_amount'] = p_con.get_ordered_amount(self.examples["project"])
                s_param['customer_reference'] = s_param.get('customer_reference', project['customer_reference'])
                s_param['uppnex_project_id'] = s_param.get('uppnex_project_id', project['uppnex_id'])
            s_param['success'] = sequencing_success(s_param, cutoffs)
            s_param.update({k:"N/A" for k in s_param.keys() if s_param[k] is None})
            make_note("{}.pdf".format(s["barcode_name"]), headers, paragraphs, **s_param)

