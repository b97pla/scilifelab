import os
import csv
import yaml
import couchdb
from couchdb.design import ViewDefinition
import unittest
import time
import logbook
import socket

from classes import PmFullTest

from ..classes import has_couchdb_installation

from scilifelab.db.statusdb import SampleRunMetricsConnection, VIEWS, flowcell_run_metrics, sample_run_metrics, project_summary, ProjectSummaryConnection, update_fn
from scilifelab.bcbio.qc import FlowcellRunMetricsParser, SampleRunMetricsParser
from scilifelab.pm.bcbio.utils import fc_id, fc_parts, fc_fullname

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
        prj = project_summary(**p)
        p_con.save(prj, key="project_id")

    #
    # def tearDownModule():
    #     db = couchdb.Server()
    #     for x in DATABASES:
    #         LOG.info("Deleting database {}".format(x))
    #         del db[x]
    

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestCouchDB(unittest.TestCase):

    def test_dbcon(self):
        s_con = SampleRunMetricsConnection(dbname="samples-test", username="u", password="p")
        samples = [s_con.get_entry(x) for x in s_con.name_view]
        for s in samples:
            print s

    def test_srm_upload(self):
        """Test upload of Sample Run Metrics"""
        # View results at http://localhost:5984/_utils/database.html?samples-test
        for fc in flowcells:
            _save_samples(os.path.join(flowcell_dir, fc))

    def test_fc_upload(self):
        for fc in flowcells:
            _save_flowcell(os.path.join(flowcell_dir, fc))

@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestQCUpload(PmFullTest):
    def test_qc_upload(self):
        """Test running qc upload to server"""
        self.app = self.make_app(argv = ['qc', 'upload-qc', flowcells[0], '--mtime',  '100'], extensions=['scilifelab.pm.ext.ext_qc'])
        self._run_app()

