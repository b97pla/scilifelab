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

DEMUX_STATS="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html xmlns:casava="http://www.illumina.com/casava/alignment" xmlns:str="http://exslt.org/strings">
<link rel="stylesheet" href="css/Reports.css" type="text/css">
<body>
<h1>Flowcell: AC003CCCXX</h1>
<h2>Barcode lane statistics</h2>
<div ID="ScrollableTableHeaderDiv"><table width="100%">
<col width="4%">
<col width="5%">
<col width="19%">
<col width="8%">
<col width="7%">
<col width="5%">
<col width="12%">
<col width="7%">
<col width="4%">
<col width="5%">
<col width="4%">
<col width="5%">
<col width="6%">
<col width="5%">
<col>
<tr>
<th>Lane</th>
<th>Sample ID</th>
<th>Sample Ref</th>
<th>Index</th>
<th>Description</th>
<th>Control</th>
<th>Project</th>
<th>Yield (Mbases)</th>
<th>% PF</th>
<th># Reads</th>
<th>% of raw clusters per lane</th>
<th>% Perfect Index Reads</th>
<th>% One Mismatch Reads (Index)</th>
<th>% of &gt;= Q30 Bases (PF)</th>
<th>Mean Quality Score (PF)</th>
</tr>
</table></div>
<div ID="ScrollableTableBodyDiv"><table width="100%">
<col width="4%">
<col width="5%">
<col width="19%">
<col width="8%">
<col width="7%">
<col width="5%">
<col width="12%">
<col width="7%">
<col width="4%">
<col width="5%">
<col width="4%">
<col width="5%">
<col width="6%">
<col width="5%">
<col>
<tr>
<td>1</td>
<td>P001_101_index3</td>
<td>hg19</td>
<td>CAGATC</td>
<td>J__Doe_00_01</td>
<td>TGACCA</td>
<td>J__Doe_00_01</td>
<td>3,942</td>
<td>100.00</td>
<td>39,034,396</td>
<td>7.94</td>
<td>92.57</td>
<td>7.43</td>
<td>90.05</td>
<td>35.22</td>
</tr>
<tr>
<td>1</td>
<td>lane1</td>
<td>unknown</td>
<td>Undetermined</td>
<td>Clusters with unmatched barcodes for lane 1</td>
<td>N</td>
<td>Undetermined_indices</td>
<td>7,867</td>
<td>100.00</td>
<td>77,892,454</td>
<td>15.85</td>
<td>0.00</td>
<td>0.00</td>
<td>87.28</td>
<td>34.16</td>
</tr>
<tr>
<td>2</td>
<td>P001_102_index6</td>
<td>hg19</td>
<td>ACAGTG</td>
<td>J__Doe_00_01</td>
<td>N</td>
<td>J__Doe_00_01</td>
<td>3,095</td>
<td>100.00</td>
<td>30,641,418</td>
<td>6.23</td>
<td>91.36</td>
<td>8.64</td>
<td>90.94</td>
<td>35.54</td>
</tr>
<tr>
<td>2</td>
<td>lane2</td>
<td>unknown</td>
<td>Undetermined</td>
<td>Clusters with unmatched barcodes for lane 2</td>
<td>N</td>
<td>Undetermined_indices</td>
<td>7,867</td>
<td>100.00</td>
<td>77,892,454</td>
<td>15.85</td>
<td>0.00</td>
<td>0.00</td>
<td>87.28</td>
<td>34.16</td>
</tr>
</table></div>
<p></p>
<h2>Sample information</h2>
<div ID="ScrollableTableHeaderDiv"><table width="100%">
<col width="10%">
<col width="10%">
<col width="7%">
<col>
<tr>
<th>Sample<p></p>ID</th>
<th>Recipe</th>
<th>Operator</th>
<th>Directory</th>
</tr>
</table></div>
<div ID="ScrollableTableBodyDiv"><table width="100%">
<col width="10%">
<col width="10%">
<col width="7%">
<col>
<tr>
<td>P364_101_index7</td>
<td>R1</td>
<td>NN</td>
<td>/srv/illumina/120924_SN0002_0003_AC003CCCXX/Unaligned/Project_D__Klevebring_12_03/Sample_P364_101_index7</td>
</tr>
<tr>
<td>P364_102_index12</td>
<td>R1</td>
<td>NN</td>
<td>/srv/illumina/120924_SN0002_0003_AC003CCCXX/Unaligned/Project_D__Klevebring_12_03/Sample_P364_102_index12</td>
</tr>
</table></div>
<p>bcl2fastq-1.8.3</p>
</body>
</html>
"""


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
        print samples_d
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
        self.app = self.make_app(argv = ['qc', 'upload-qc', flowcells[0], '--mtime',  '10000'], extensions=['scilifelab.pm.ext.ext_qc', 'scilifelab.pm.ext.ext_couchdb'])
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
        self.assertEqual(fc["samplesheet_csv"][2]["SampleID"], "P001_101_index3")        

    def test_qc_upload(self):
        """Test running qc upload to server"""
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
