import os
import yaml
import couchdb
from couchdb.design import ViewDefinition
import unittest
import logbook
import xml.etree.cElementTree as ET
import shutil

from classes import PmFullTest
from ..classes import has_couchdb_installation

from scilifelab.db.statusdb import SampleRunMetricsConnection, VIEWS, ProjectSummaryDocument, ProjectSummaryConnection, FlowcellRunMetricsConnection
from scilifelab.bcbio.qc import FlowcellRunMetricsParser, SampleRunMetricsParser,  XmlToDict

filedir = os.path.dirname(os.path.abspath(__file__))
dirs = {'production': os.path.join(filedir, "data", "production")}

LOG = logbook.Logger(__name__)

flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
flowcell_dir = os.path.join(filedir, "data", "archive")
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

def setUpModule():
    """Create test databases in local server"""
    if not has_couchdb:
        return
    server = couchdb.Server()

    ## Create databases
    for x in DATABASES:
        if not server.__contains__(x):
            LOG.info("Creating database {}".format(x))
            server.create(x)
    ## Create views for flowcells and samples
    for dbname in DATABASES:
        dblab = dbname.replace("-test", "")
        db = server[dbname]
        for k,v in VIEWS[dblab].items():
            for title, view in v.items():
                viewdef = ViewDefinition(k, title, view)
                viewdef.sync(db)
    
    ## Create and upload project summary
    with open(os.path.join(filedir, "data", "config", "project_summary.yaml")) as fh:
        prj_sum = yaml.load(fh)
    db = server["samples-test"]
    p_con = ProjectSummaryConnection(dbname="projects-test", username="u", password="p")
    for p in prj_sum:
        prj = ProjectSummaryDocument(**p)
        p_con.save(prj, key="project_name")


    #
    # def tearDownModule():
    #     db = couchdb.Server()
    #     for x in DATABASES:
    #         LOG.info("Deleting database {}".format(x))
    #         del db[x]
    

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestCouchDB(unittest.TestCase):

    def test_dbcon(self):
        """Test database connection and that we get expected values."""
        s_con = SampleRunMetricsConnection(dbname="samples-test", username="u", password="p")
        samples = [s_con.get_entry(x) for x in s_con.name_view]
        samples_d = {x["name"]:x for x in samples}
        self.assertEqual(samples_d["1_120924_AC003CCCXX_TGACCA"]["date"], "120924")
        self.assertEqual(samples_d["1_121015_BB002BBBXX_TGACCA"]["flowcell"], "BB002BBBXX")
        self.assertEqual(samples_d["2_120924_AC003CCCXX_ACAGTG"]["entity_type"], "sample_run_metrics")
        self.assertEqual(samples_d["3_120924_AC003CCCXX_ACAGTG"]["lane"], "3")
        self.assertEqual(samples_d["4_120924_AC003CCCXX_CGTTAA"]["sequence"], "CGTTAA")
        self.assertEqual(samples_d["2_121015_BB002BBBXX_TGACCA"]["project_id"], "P002")
        fc_con = FlowcellRunMetricsConnection(dbname="flowcells-test", username="u", password="p")
        flowcells = [fc_con.get_entry(x) for x in fc_con.name_view]
        flowcells_d = {x["name"]:x for x in flowcells}
        self.assertEqual(flowcells_d["120924_AC003CCCXX"]["name"], "120924_AC003CCCXX")
        self.assertEqual(flowcells_d["121015_BB002BBBXX"]["name"], "121015_BB002BBBXX")
        self.assertEqual(flowcells_d["120924_AC003CCCXX"]["entity_type"], "flowcell_run_metrics")
        p_con = ProjectSummaryConnection(dbname="projects-test", username="u", password="p")
        projects = [p_con.get_entry(x) for x in p_con.name_view]
        projects_d = {x["project_name"]:x for x in projects}
        self.assertEqual(projects_d["J.Doe_00_01"]["min_m_reads_per_sample_ordered"], 0.1)
        self.assertEqual(projects_d["J.Doe_00_01"]["no_of_samples"], 2)        
        self.assertEqual(set(projects_d["J.Doe_00_01"]["samples"].keys()),set(["P001_101_index3","P001_102","P001_103"]))
        self.assertEqual(projects_d["J.Doe_00_01"]["customer_reference"], "GnuGenome")
        self.assertEqual(projects_d["J.Doe_00_02"]["min_m_reads_per_sample_ordered"], 0.2)        
        self.assertEqual(projects_d["J.Doe_00_03"]["samples"].keys(),["3_index6"])
        self.assertIn("A", projects_d["J.Doe_00_03"]["samples"]["3_index6"]["library_prep"])

 

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestQCUpload(PmFullTest):
    def setUp(self):
        """FIXME: All other tests depend on data being uploaded, so
        these are not real unit tests. The setup to TestQCUpload has to
        be run prior to other tests, else unexpected failures will
        occur."""
        self.app = self.make_app(argv = ['qc', 'upload-qc', flowcells[0], '--mtime',  '10000'], extensions=['scilifelab.pm.ext.ext_qc', 'scilifelab.pm.ext.ext_couchdb'])
        self._run_app()
        self.app = self.make_app(argv = ['qc', 'upload-qc', flowcells[1], '--mtime',  '10000'], extensions=['scilifelab.pm.ext.ext_qc', 'scilifelab.pm.ext.ext_couchdb'])
        self._run_app()
        self.s_con = SampleRunMetricsConnection(dbname="samples-test", username="u", password="p")
        self.p_con = ProjectSummaryConnection(dbname="projects-test", username="u", password="p")
        self.fc_con = FlowcellRunMetricsConnection(dbname="flowcells-test", username="u", password="p")

    def test_samplesheet(self):
        """Test samplesheet upload"""
        fc = self.fc_con.get_entry("120924_AC003CCCXX")
        self.assertEqual(fc["samplesheet_csv"][0]["Index"], "TGACCA")
        self.assertEqual(fc["samplesheet_csv"][0]["Description"], "J__Doe_00_01")
        self.assertEqual(fc["samplesheet_csv"][0]["FCID"], "C003CCCXX")
        self.assertEqual(fc["samplesheet_csv"][1]["SampleRef"], "hg19")
        self.assertEqual(fc["samplesheet_csv"][2]["SampleID"], "P002_101_index3")

    def test_qc_upload(self):
        """Test running qc upload to server. Slightly circular testing
        here - I setup the module with qc update so by definition the
        test must 'work'"""
        self.app = self.make_app(argv = ['qc', 'upload-qc', flowcells[1], '--mtime',  '100'], extensions=['scilifelab.pm.ext.ext_qc',  'scilifelab.pm.ext.ext_couchdb'])
        self._run_app()
        s = self.s_con.get_entry("4_120924_AC003CCCXX_CGTTAA")
        self.assertIsNone(s["project_sample_name"])
        self.assertEqual(s["project_id"], "P003")
        
    def test_qc_update(self):
        """Test running qc update of a project id"""
        s = self.s_con.get_entry("4_120924_AC003CCCXX_CGTTAA")
        s["project_id"]= None
        self.assertIsNone(s["project_id"])
        self.s_con.save(s)
        self.app = self.make_app(argv = ['qc', 'update', '--sample_prj', projects[2], '--project_id', 'P003', '--debug', '--force'], extensions=['scilifelab.pm.ext.ext_qc', 'scilifelab.pm.ext.ext_couchdb'])
        self._run_app()
        s = self.s_con.get_entry("4_120924_AC003CCCXX_CGTTAA")
        self.assertEqual(s["project_id"], "P003")

    def test_qc_update_sample_names(self):
        """Test running qc update of project sample names"""
        s1 = self.s_con.get_entry("1_120924_AC003CCCXX_TGACCA")
        s2 = self.s_con.get_entry("2_120924_AC003CCCXX_ACAGTG")
        s1["project_sample_name"] = None
        s2["project_sample_name"] = None
        self.assertIsNone(s1["project_sample_name"])
        self.assertIsNone(s2["project_sample_name"])
        self.s_con.save(s1)
        self.s_con.save(s2)
        sample_map = {'P001_101_index3': 'P001_101_index3', 'P001_102_index6':'P001_102'}
        self.app = self.make_app(argv = ['qc', 'update', '--sample_prj', projects[0], '--names', "{}".format(sample_map), '--debug', '--force'], extensions=['scilifelab.pm.ext.ext_qc',  'scilifelab.pm.ext.ext_couchdb'])
        self._run_app()
        s1 = self.s_con.get_entry("1_120924_AC003CCCXX_TGACCA")
        s2 = self.s_con.get_entry("2_120924_AC003CCCXX_ACAGTG")
        self.assertEqual(s1["project_sample_name"], "P001_101_index3")
        self.assertEqual(s2["project_sample_name"], "P001_102")

class TestMetricsParser(PmFullTest):
    def setUp(self):
        self.sample_kw = dict(flowcell="AC003CCCXX", date="120924", lane=1, barcode_name='P001_101_index3', sample_prj="J.Doe_00_01".replace("__", "."), barcode_id="1", sequence="TGACCA")
        self.fc_kw = dict(fc_date = "120924", fc_name = "AC003CCCXX")
        self.fcdir = os.path.join(flowcell_dir, "120924_SN0002_0003_AC003CCCXX")
        # Add Demultiplex_stats for j_doe_00_01
        demux_stats_file = os.path.join(filedir, "data", "archive",  "120924_SN0002_0003_AC003CCCXX", "Unaligned", "Basecall_Stats_C003CCCXX", "Demultiplex_Stats.htm")
        if not os.path.exists(demux_stats_file):
            # TODO: use links...
            with open(os.path.join("data", "db", "demux_stats.htm")) as fh:
                DEMUX_STATS = fh.readlines()
            with open(demux_stats_file, "w") as fh:
                fh.write(DEMUX_STATS)

    def test_get_bc_count(self):
        parser = SampleRunMetricsParser(os.path.join(project_dir, "J.Doe_00_01", "P001_101_index3", "120924_AC003CCCXX"))
        bc_count = parser.get_bc_count(**self.sample_kw)
        self.assertEqual(bc_count, 0)

    def test_get_bc_count_demux_stats(self):
        parser = SampleRunMetricsParser(os.path.join(project_dir, "J.Doe_00_01", "P001_101_index3", "120924_AC003CCCXX"))
        bc_count = parser.get_bc_count(**self.sample_kw)
        fc_parser = FlowcellRunMetricsParser(self.fcdir)
        data = fc_parser.parse_demultiplex_stats_htm(**self.fc_kw)
        bc_count = parser.get_bc_count(demultiplex_stats=data, **self.sample_kw)
        self.assertEqual(str(bc_count), str(19517198))

    def test_parseRunParameters(self):
        parser = FlowcellRunMetricsParser(self.fcdir)
        data = parser.parseRunParameters(**self.fc_kw)
        self.assertEqual(data['Setup']['FPGADynamicFocusSettings']['CVGainPosLocked'],'500')
        self.assertEqual(data['Setup']['Reads']['Read'][0]['NumCycles'],'101')
        self.assertEqual(data['Setup']['SelectedSections']['Section'][1]['Name'],'B_1')

    def test_parse_demux_stats(self):
        """Test parsing of a Demultiplex_Stats.htm file"""
        parser = FlowcellRunMetricsParser(self.fcdir)
        data = parser.parse_demultiplex_stats_htm(**self.fc_kw)
        self.assertEqual(data['Barcode_lane_statistics'][0]['# Reads'], '39,034,396')

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestDbConnection(unittest.TestCase):
    def setUp(self):
        self.user = "user"
        self.pw = "pw"
        self.url = "localhost"
        self.examples = {"sample":"1_120924_AC003CCCXX_TGACCA",
                         "flowcell":"AC003CCCXX",
                         "project":"J.Doe_00_01"}
        self.p_con = ProjectSummaryConnection(dbname="projects-test", username=self.user, password=self.pw, url=self.url)

    def test_connection(self):
        """Test database connection"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        self.assertEqual(sample_con.url_string, "http://{}:5984".format(self.url))

    def test_get_flowcell(self):
        """Test getting a flowcell for a given sample"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        fc = sample_con.get_entry(self.examples["sample"], "flowcell")
        self.assertEqual(str(fc), self.examples["flowcell"])

    def test_get_sample_ids(self):
        """Test getting sample ids given flowcell and sample_prj"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"])
        LOG.info("Number of samples before subsetting: " + str(len(sample_ids)))
        self.assertEqual(len(sample_ids), 5)
        sample_ids = sample_con.get_sample_ids(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info( "Number of samples after subsetting: " + str(len(sample_ids)))
        self.assertEqual(len(sample_ids), 2)

    def test_get_samples(self):
        """Test getting samples given flowcell and sample_prj."""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)

        samples = sample_con.get_samples(fc_id=self.examples["flowcell"])
        LOG.info("Selecting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 5)
        samples = sample_con.get_samples(fc_id=self.examples["flowcell"], sample_prj=self.examples["project"])
        LOG.info("Selecting on flowcell, subsetting on project: " + str(len(samples)))
        self.assertEqual(len(samples), 2)

        samples = sample_con.get_samples(sample_prj=self.examples["project"])
        LOG.info("Selecting on project: " + str(len(samples)))
        self.assertEqual(len(samples), 3)
        samples = sample_con.get_samples(sample_prj=self.examples["project"], fc_id=self.examples["flowcell"])
        LOG.info("Selecting on project, subsetting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 2)

    def test_get_samples_wrong_info(self):
        """Test getting samples when either flowcell or project id information is wrong"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)

        samples = sample_con.get_samples(sample_prj="bogusproject", fc_id=self.examples["flowcell"])
        LOG.info("Selecting on bogus project, subsetting on flowcell: " + str(len(samples)))
        self.assertEqual(len(samples), 0)
        
                
    def test_get_project_sample_ids(self):
        """Test getting project sample ids"""
        sample_con = SampleRunMetricsConnection(dbname="samples-test", username=self.user, password=self.pw, url=self.url)
        sample_ids = sample_con.get_sample_ids(sample_prj=self.examples["project"])
        sample_names = [sample_con.db.get(x)["name"] for x in sample_ids]
        self.assertEqual(set(sample_names) , set(['1_120924_AC003CCCXX_TGACCA', '2_120924_AC003CCCXX_ACAGTG', '1_121015_BB002BBBXX_TGACCA']))
        
    def test_get_latest_library_prep(self):
        """Test getting latest library prep"""
        prj = self.p_con.get_entry("J.Doe_00_01")
        prj['samples']['P001_102']['library_prep']['B'] = {'sample_run_metrics': {'2_120924_AC003CCCXX_TTGGAA': None}}
        self.p_con.save(prj)
        preps = self.p_con.get_latest_library_prep(project_name=self.examples["project"])
        srm = [x for l in preps.values() for x in l]
        # Make sure A prep not in list
        self.assertNotIn('2_120924_AC003CCCXX_ACAGTG', srm)
        # Make sure B prep in list
        self.assertIn('2_120924_AC003CCCXX_TTGGAA', srm)
        # Reset data
        prj = self.p_con.get_entry("J.Doe_00_01")
        del prj['samples']['P001_102']['library_prep']['B']
        self.p_con.save(prj)

    def test_get_barcode_lane_statistics(self):
        """Test getting barcode lane statistics from flowcell database"""
        fc_con = FlowcellRunMetricsConnection(dbname="flowcells-test", username="u", password="p")
        # Try getting wrong sample name, should return None
        data = fc_con.get_barcode_lane_statistics("J.Doe_00_01", "P001_101_index6", "120924_AC003CCCXX", "1")
        self.assertEqual(data, (None, None))
        data = fc_con.get_barcode_lane_statistics("J.Doe_00_01", "P001_101_index3", "120924_AC003CCCXX", "1")
        self.assertEqual(data, (u'35.22', u'90.05'))
