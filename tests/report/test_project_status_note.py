import os
import re
import unittest
import itertools
import ConfigParser
import logbook
from scilifelab.report.rl import make_example_project_note, make_note, project_note_paragraphs, project_note_headers, make_sample_table
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectSummaryConnection

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
LOG = logbook.Logger(__name__)


## Cutoffs
cutoffs = {
    "phix_err_cutoff" : 2.0,
    }

## parameters
parameters = {
    "project_name" : None,
    "customer_reference": None,
    "uppnex_project_id" : None,
    "finished":"Not finished, or cannot yet assess if finished.",
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
            LOG.warning("No such file {}; will not run database connection tests".format(os.path.join(os.getenv("HOME"), "dbcon.ini")))
        else:
            config = ConfigParser.ConfigParser()
            config.readfp(open(os.path.join(os.getenv("HOME"), "dbcon.ini")))
            self.url = config.get("couchdb", "url")
            self.user = config.get("couchdb", "username")
            self.pw = config.get("couchdb", "password")
            self.examples = {"sample":config.get("examples", "sample"),
                             "flowcell":config.get("examples", "flowcell"),
                             "project":config.get("examples", "project")}

    def test_make_example_project_note(self):
        """Make example project note"""
        make_example_project_note(os.path.join(filedir, "project_test.pdf"))

    def test_make_project_note(self):
        """Make a project note subset by flowcell and project"""
        if not self.examples:
            LOG.info("Not running test")
            return
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        fc_con = FlowcellRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        paragraphs = project_note_paragraphs()
        headers = project_note_headers()
        param = parameters
        project = p_con.get_entry(self.examples["project"])
        if not project:
            print "No project named {}".format(self.examples["project"])
            return
        if project:
            ordered_amount = p_con.get_ordered_amount(self.examples["project"])
        else:
            return
            ordered_amount = self.pargs.ordered_million_reads
            

        ## Start collecting the data
        sample_table = []
        sample_list = project['samples']
        param.update({key:project.get(ps_to_parameter[key], None) for key in ps_to_parameter.keys()})
        samples = p_con.map_srm_to_name(self.examples["project"], check_consistency=True)
        all_passed = True
        for k,v in samples.items():
            print k, v
            if re.search("Unexpected", k):
                continue
            project_sample = sample_list[v['sample']]
            vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
            vals['Status'] = project_sample.get("status", "N/A")
            vals['MOrdered'] = ordered_amount
            vals['BarcodeSeq'] = s_con.get_entry(k, "sequence")
            vals.update({k:"N/A" for k in vals.keys() if vals[k] is None})
            if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
            sample_table.append([vals[k] for k in table_keys])
        if all_passed: param["finished"] = 'Project finished.'
        sample_table.sort()
        sample_table = list(sample_table for sample_table,_ in itertools.groupby(sample_table))
        sample_table.insert(0, ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
        paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
        make_note("{}.pdf".format(self.examples["project"]), headers, paragraphs, **param)

    def test_sample_map(self):
        """Test getting a sample mapping"""
        if not self.examples:
            LOG.info("Not running test")
            return
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        sample_map = p_con.map_name_to_srm(self.examples["project"], use_ps_map=False, check_consistency=True)
        print sample_map

    def test_sample_map_fc(self):
        """Test getting a sample mapping subset by flowcell"""
        if not self.examples:
            LOG.info("Not running test")
            return
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        sample_map = p_con.map_name_to_srm(self.examples["project"], self.examples["flowcell"], use_ps_map=False, check_consistency=True)
        print sample_map

    def test_srm_map(self):
        """Test getting a sample mapping from srm to project samples"""
        if not self.examples:
            LOG.info("Not running test")
            return
        p_con = ProjectSummaryConnection(username=self.user, password=self.pw, url=self.url)
        samples = p_con.map_srm_to_name(self.examples["project"], fc_id=self.examples["flowcell"], use_ps_map=False, check_consistency=True)
        print samples

    def test_srm_sample_map(self):
        """Make sample map from sample run metrics"""
        if not self.examples:
            LOG.info("Not running test")
            return
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = s_con.get_samples(sample_prj=self.examples["project"])
        sample_names = {x["name"]:x["_id"] for x in samples}
        print len(sample_names.keys())
        print len(list(set(sample_names.keys())))
        print sample_names

    def test_name_view(self):
        """Test unique ids in name view"""
        if not self.examples:
            LOG.info("Not running test")
            return
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        name_view  = s_con.db.view("names/name", reduce=False)
        print s_con.name_view.keys()[0:10]
        print len(s_con.name_view.keys())
        print len(list(set(s_con.name_view.keys())))
        print len(name_view)

