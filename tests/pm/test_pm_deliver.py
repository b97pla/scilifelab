"""
Test pm-delivery application
"""
import os
import ast
import logbook
import unittest
import copy
from cement.core import handler
from test_default import PmTest
from scilifelab.pm.core.deliver import *
from scilifelab.db.statusdb import SampleRunMetricsConnection, SampleRunMetricsDocument

from ..classes import has_couchdb_installation

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class PmProductionTest(PmTest):
    def setUp(self):
        self.url = "localhost"
        self.user = "user"
        self.pw = "pw"
        # ['1_120924_AC003CCCXX_TGACCA', '2_120924_AC003CCCXX_ACAGTG', '1_121015_BB002BBBXX_TGACCA']
        self.examples = {"sample": "1_120924_AC003CCCXX_TGACCA",
                         "flowcell": "AC003CCCXX",
                         "project":"J.Doe_00_01", 
                         "exclude_samples": "{'P001_102':[]}"}

    @classmethod
    @unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
    def setUpClass(cls):
        # Temporarily add new sample for use in exclusion tests
        s_con = SampleRunMetricsConnection(dbname="samples-test", username="u", password="p")
        s = s_con.get_entry("1_121015_BB002BBBXX_TGACCA")
        kw = copy.deepcopy(s)
        del kw["_id"]
        new_s = SampleRunMetricsDocument(**kw)
        new_s["sequence"] = "AGTTGA"
        new_s["name"] = "1_121015_BB002BBBXX_AGTTGA"
        s_con.save(new_s)

        kw = copy.deepcopy(s)
        del kw["_id"]
        new_s = SampleRunMetricsDocument(**kw)
        new_s["sample_prj"] = "j-doe_00_01"
        new_s["sequence"] = "CGAACG"
        new_s["name"] = "1_121015_BB002BBBXX_CGAACG"
        s_con.save(new_s)

        s = s_con.get_entry("3_120924_AC003CCCXX_ACAGTG")
        kw = copy.deepcopy(s)
        del kw["_id"]
        new_s = SampleRunMetricsDocument(**kw)
        new_s["sample_prj"] = "j-doe_00_02"
        new_s["sequence"] = "GGAAGG"
        new_s["name"] = "3_120924_AC003CCCXX_GGAAGG"
        s_con.save(new_s)

    @classmethod
    @unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
    def tearDownClass(cls):
        s_con = SampleRunMetricsConnection(dbname="samples-test", username="u", password="p")
        s = s_con.get_entry("1_121015_BB002BBBXX_AGTTGA")
        doc = s_con.db.get(s["_id"])
        s_con.db.delete(doc)
        s = s_con.get_entry("1_121015_BB002BBBXX_CGAACG")
        doc = s_con.db.get(s["_id"])
        s_con.db.delete(doc)
        s = s_con.get_entry("3_120924_AC003CCCXX_GGAAGG")
        doc = s_con.db.get(s["_id"])
        s_con.db.delete(doc)

    # Will currently fail since no PhiX in document
    def test_sample_status(self):
        self.app = self.make_app(argv = ['report', 'sample-status', self.examples["project"], self.examples["flowcell"], '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(data['s_param']['P001_101_index3']['scilifelab_name'], 'P001_101_index3')
        self.assertEqual(data['s_param']['P001_101_index3']['customer_reference'], 'GnuGenome')
        self.assertEqual(len(data['sample_runs'].keys()), 2)

    def test_project_status(self):
        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"]],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(data['param']['project_name'], 'J.Doe_00_01')
        self.assertEqual(data['param']['ordered_amount'], 0.1)
        
    def test_sample_status_custom(self):
        """Test sample status note generation with command line customizations"""
        self.app = self.make_app(argv = ['report', 'sample_status', self.examples["project"], self.examples["flowcell"], '--debug', '--customer_reference', 'MyCustomerReference', '--uppnex_id', 'MyUppnexID', '--ordered_million_reads', '10'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(data['s_param']['P001_101_index3']['uppnex_project_id'], 'MyUppnexID')
        self.assertEqual(data['s_param']['P001_101_index3']['customer_reference'], 'MyCustomerReference')
        self.assertEqual(data['s_param']['P001_101_index3']['ordered_amount'], 10)

    def test_project_status_exclude_samples(self):
        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"], '--debug',  '--exclude_sample_ids', "{'P001_102':[]}"],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        samples = [x[0] for x in  data['table']]
        self.assertNotIn("P001_102", samples)

        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"], '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"], '--debug',  '--exclude_sample_ids', "{'P001_101_index3':['AGTTGA']}"],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        barcodes = [x[2] for x in  data['table']]
        self.assertNotIn("AGTTGA", barcodes)

    def test_ordered_amount(self):
        """Test setting ordered amount to different values for different samples"""
        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"], '--debug',  '-o', "{'P001_101_index3':10, 'P001_102':20}"],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        ordered = {x[0]:x[4] for x in data['table']}
        self.assertEqual(ordered["P001_101_index3"], 10)
        self.assertEqual(ordered["P001_102"], 20)

    def test_bc_count(self):
        """Test setting ordered amount to different values for different samples"""
        self.app = self.make_app(argv = ['report', 'sample_status', self.examples["project"], self.examples["flowcell"],  '--debug',  '--phix', '0.1', '--bc_count', "{'P001_101_index3':1300000, 'P001_102_index6':20000}"],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(data['s_param']['P001_101_index3']['ordered_amount'], 0.1)
        self.assertEqual(data['s_param']['P001_102_index6']['ordered_amount'], 0.1)
        self.assertEqual(data['s_param']['P001_101_index3']['rounded_read_count'], 1.3)
        self.assertEqual(data['s_param']['P001_102_index6']['rounded_read_count'], 0.02)
        self.assertEqual(data['s_param']['P001_101_index3']['success'], 'Successful run.')
        self.assertEqual(data['s_param']['P001_102_index6']['success'], 'The yield may be lower than expected.')

    def test_sample_aliases(self):
        """Test setting sample aliases to different values for different samples"""
        self.app = self.make_app(argv = ['report', 'project_status', 'J.Doe_00_03', '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        # This should fail since P003_101_index6 != 3_index6
        self.assertEqual(len(data['table']), 2)

        self.app = self.make_app(argv = ['report', 'project_status', 'J.Doe_00_03', '--sample_alias', "{'P003_101_index6':'3_index6'}", '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        samples = [x[0] for x in  data['table']]
        self.assertIn("3_index6", samples)
        
    def test_project_alias(self):
        """Test setting project alias"""
        self.app = self.make_app(argv = ['report', 'project_status', 'J.Doe_00_01', '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        barcodes = [x[2] for x in data['table'][1:]]
        self.assertEqual(len(data['table']), 5)
        self.assertEqual(set(['AGTTGA', 'TGACCA', 'ACAGTG', 'N/A']), set(barcodes))
        
        self.app = self.make_app(argv = ['report', 'project_status', 'J.Doe_00_01', '--project_alias', '["j-doe_00_01"]', '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        barcodes = [x[2] for x in data['table'][1:]]
        self.assertEqual(len(data['table']), 6)
        self.assertIn("CGAACG", barcodes)

        self.app = self.make_app(argv = ['report', 'sample-status', 'J.Doe_00_02', 'AC003CCCXX', '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(len(data['sample_runs'].keys()), 2)
        self.assertNotIn("3_120924_AC003CCCXX_GGAAGG", data['sample_runs'].keys())

        self.app = self.make_app(argv = ['report', 'sample-status', 'J.Doe_00_02', 'AC003CCCXX', '--project_alias', '["j-doe_00_02"]', '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
        data = ast.literal_eval(self.app._output_data['debug'].getvalue())
        self.assertEqual(len(data['sample_runs'].keys()), 3)
        self.assertIn("3_120924_AC003CCCXX_GGAAGG", data['sample_runs'].keys())

