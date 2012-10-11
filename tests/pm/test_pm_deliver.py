"""
Test pm-delivery application
"""
import os
import ConfigParser
from cement.core import handler
from test_default import PmTest
from scilifelab.pm.core.deliver import *

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

class PmProductionTest(PmTest):
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

        
    def test_1_sample_status(self):
        self.app = self.make_app(argv = ['report', 'sample_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], self.examples["flowcell"], '--debug', '--use_bc_map'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_1_project_status(self):
        self.app = self.make_app(argv = ['report', 'project_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], '--debug', '--use_bc_map'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_3_sample_status_ps_map(self):
        self.app = self.make_app(argv = ['report', 'sample_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], self.examples["flowcell"], '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_3_project_status_ps_map(self):
        self.app = self.make_app(argv = ['report', 'project_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], '--debug'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()


    def test_2_sample_status_custom(self):
        """Test sample status note generation with command line customizations"""
        self.app = self.make_app(argv = ['report', 'sample_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], self.examples["flowcell"], '--debug', '--customer_reference', 'MyCustomerReference', '--uppnex_id', 'MyUppnexID', '--ordered_million_reads', '10', '--use_bc_map'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()

    def test_2_project_status_custom(self):
        """Test project status note generation with command line customizations"""
        self.app = self.make_app(argv = ['report', 'project_status', '--user', self.user, '--password', self.pw, '--url', self.url, self.examples["project"], '--debug', '--customer_reference', 'MyCustomerReference', '--uppnex_id', 'MyUppnexID', '--use_bc_map', '--ordered_million_reads', '1'],extensions=['scilifelab.pm.ext.ext_couchdb'])
        handler.register(DeliveryReportController)
        self._run_app()
