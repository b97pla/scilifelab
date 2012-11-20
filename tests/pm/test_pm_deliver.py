"""
Test pm-delivery application
"""
import os
import logbook
from cement.core import handler
from test_default import PmTest
from scilifelab.pm.core.deliver import *

from ..classes import has_couchdb_installation

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

class PmProductionTest(PmTest):
    def setUp(self):
        self.url = "localhost"
        self.user = "user"
        self.pw = "pw"
        # ['1_120924_AC003CCCXX_TGACCA', '2_120924_AC003CCCXX_ACAGTG', '1_121015_BB002BBBXX_TGACCA']
        self.examples = {"sample": "1_120924_AC003CCCXX_TGACCA",
                         "flowcell": "AC003CCCXX",
                         "project":"J.Doe_00_01", 
                         "exclude_samples": "{'P001_102_index6':[]}"}

    # Will currently fail since no PhiX in document
    def test_sample_status(self):
        self.app = self.make_app(argv = ['report', 'sample-status', self.examples["project"], self.examples["flowcell"], '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_project_status(self):
        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"]],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_sample_status_custom(self):
        """Test sample status note generation with command line customizations"""
        self.app = self.make_app(argv = ['report', 'sample_status', self.examples["project"], self.examples["flowcell"], '--debug', '--customer_reference', 'MyCustomerReference', '--uppnex_id', 'MyUppnexID', '--ordered_million_reads', '10'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_project_status_exclude_samples(self):
        self.app = self.make_app(argv = ['report', 'project_status', self.examples["project"], '--debug',  '--exclude_sample_ids', "{}".format(self.examples["exclude_samples"])],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
