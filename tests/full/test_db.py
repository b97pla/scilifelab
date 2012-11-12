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
from scilifelab.pm.ext.ext_couchdb import CouchdbCommandHandler
from scilifelab.pm.ext.ext_qc import update_fn
from scilifelab.db.statusdb import SampleRunMetricsConnection
from scilifelab.bcbio.qc import FlowcellRunMetrics, SampleRunMetrics
from scilifelab.pm.bcbio.utils import validate_fc_directory_format, fc_id, fc_parts, fc_fullname

filedir = os.path.dirname(os.path.abspath(__file__))
dirs = {'production': os.path.join(filedir, "data", "production")}

LOG = logbook.Logger(__name__)

flowcell = os.path.join(filedir, "data", "archive", "120924_SN0002_0003_AC003CCCXX")

project = os.path.join(filedir, "data", "production", "J.Doe_00_01")

## Try connecting to server
has_couchdb = True
try:
    server = couchdb.Server()
    dbstats = server.stats()
except socket.error as e:
    has_couchdb = False
    LOG.info("To run db tests setup a local couchdb server at http://localhost:5984")
    time.sleep(3)
    pass
    
DATABASES = ["samples-test", "projects-test", "flowcell-test"]
@unittest.skipIf(not has_couchdb, "No couchdb server running in http://localhost:5984")
class TestCouchDB(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create test databases in local server"""
        if not has_couchdb:
            return
        server = couchdb.Server()
        ## Create databases
        for x in DATABASES:
            if not server.__contains__(x):
                LOG.info("Creating database {}".format(x))
                server.create(x)
        ## Create views
        db = server["samples-test"]
        samples_name_view = ViewDefinition('names', 'name', '''function(doc) {if (!doc["name"].match(/_[0-9]+$/)) {emit(doc["name"], null);}}''')
        samples_name_view.sync(db)
        ## Create project summary
        with open(os.path.join(filedir, "data", "config", "project_summary.yaml")) as fh:
            prj_sum = yaml.load(fh)
        print prj_sum


    # @classmethod
    # def tearDownClass(cls):
    #     db = couchdb.Server()
    #     for x in DATABASES:
    #         LOG.info("Deleting database {}".format(x))
    #         del db[x]

    def test_srm_upload(self):
        """Test upload of Sample Run Metrics"""
        server = couchdb.Server()
        db = server["samples-test"]
        fcdir = flowcell
        (fc_date, fc_name) = fc_parts(flowcell)
        runinfo_csv = os.path.join(os.path.abspath(flowcell), "{}.csv".format(fc_id(flowcell)))
        try:
            with open(runinfo_csv) as fh:
                runinfo_reader = csv.reader(fh)
                runinfo = [x for x in runinfo_reader]
        except IOError as e:
            self.app.log.warn(str(e))
            raise e
        for sample in runinfo[1:]:
            d = dict(zip(runinfo[0], sample))
            sampledir = os.path.join(os.path.abspath(dirs["production"]), d['SampleProject'].replace("__", "."), d['SampleID'])
            if not os.path.exists(sampledir):
                LOG.warn("No such sample directory: {}".format(sampledir))
                continue
            sample_fcdir = os.path.join(sampledir, fc_fullname(flowcell))
            if not os.path.exists(sample_fcdir):
                LOG.warn("No such sample flowcell directory: {}".format(sample_fcdir))
                continue
            runinfo_yaml_file = os.path.join(sample_fcdir, "{}-bcbb-config.yaml".format(d['SampleID']))
            if not os.path.exists(runinfo_yaml_file):
                LOG.warn("No such yaml file for sample: {}".format(runinfo_yaml_file))
                raise IOError(2, "No such yaml file for sample: {}".format(runinfo_yaml_file), runinfo_yaml_file)
            with open(runinfo_yaml_file) as fh:
                runinfo_yaml = yaml.load(fh)
            if not runinfo_yaml['details'][0].get("multiplex", None):
                LOG.warn("No multiplex information for sample {}".format(d['SampleID']))
                continue
            sample_kw = dict(path=sample_fcdir, flowcell=fc_name, date=fc_date, lane=d['Lane'], barcode_name=d['SampleID'], sample_prj=d['SampleProject'].replace("__", "."), barcode_id=runinfo_yaml['details'][0]['multiplex'][0]['barcode_id'], sequence=runinfo_yaml['details'][0]['multiplex'][0]['sequence'])
            obj = SampleRunMetrics(**sample_kw)
            obj.read_picard_metrics()
            obj.parse_fastq_screen()
            obj.parse_bc_metrics()
            obj.read_fastqc_metrics()
            db.save(obj)
            # View results at http://localhost:5984/_utils/database.html?samples-test
            

    def test_fc_upload(self):
        server = couchdb.Server()
        db = server["flowcell-test"]
        fcdir = flowcell
        (fc_date, fc_name) = fc_parts(flowcell)
        fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
        fcobj = FlowcellRunMetrics(**fc_kw)
        fcobj.parse_illumina_metrics(fullRTA=False)
        fcobj.parse_bc_metrics()
        fcobj.parse_demultiplex_stats_htm()
        fcobj.parse_samplesheet_csv()
        db.save(fcobj)
