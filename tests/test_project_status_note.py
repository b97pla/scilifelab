import os
import unittest
import ConfigParser
from scilifelab.report import sequencing_success
from scilifelab.report.rl import make_example_project_note, make_note, project_note_paragraphs, project_note_headers, make_sample_table
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectQCSummaryConnection, ProjectSummaryConnection

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
    "finished":None,
}

## key mapping from sample_run_metrics to parameter keys
srm_to_parameter = {"project_name":"sample_prj"}
## mapping project_summary to parameter keys
ps_to_parameter = {"customer_reference":"customer_reference", "uppnex_project_id":"uppnex_id", "scilife_name":"scilife_name", "customer_name":"customer_name", "project_name":"project_id"}
## mapping project sample to table
srm_to_table = {'BarcodeSeq':'sequence'}
table_keys = ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status']
prjs_to_table = {'ScilifeID':'scilife_name', 'CustomerID':'customer_name', 'MSequenced':'m_reads_sequenced', 'MOrdered':'min_m_reads_per_sample_ordered', 'Status':'status'}

class TestProjectStatusNote(unittest.TestCase):
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

    def test_1_make_example_project_note(self):
        """Make example project note"""
        make_example_project_note(os.path.join(filedir, "project_test.pdf"))

    def test_2_make_project_note(self):
        """Make a project note subset by flowcell and project"""
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc_con = FlowcellRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        paragraphs = project_note_paragraphs()
        headers = project_note_headers()
        param = parameters
        project = p_con.get_entry(self.examples["project"])
        ## Start collecting the data
        sample_table = []
        sample_table.append(['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
        if project:
            sample_list = project['samples']
            param.update({key:project.get(ps_to_parameter[key], None) for key in ps_to_parameter.keys()})
            sample_map = p_con.map_sample_run_names(self.examples["project"], self.examples["flowcell"])
            for k,v in sample_map.items():
                project_sample = sample_list[v['project_sample']]
                vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
                vals['BarcodeSeq'] = s_con.get_entry(k, "sequence")
                vals.update({k:"N/A" for k in vals.keys() if vals[k] is None})
                sample_table.append([vals[k] for k in table_keys])
            paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
            make_note("{}.pdf".format(self.examples["project"]), headers, paragraphs, **param)
