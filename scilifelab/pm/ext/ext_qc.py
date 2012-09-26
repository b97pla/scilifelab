"""
ext_qc.py
"""
import os
import re
import csv
import yaml
import couchdb
from datetime import datetime
import time
from scilifelab.utils.timestamp import utc_time

from cement.core import backend, controller, handler, hook
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.timestamp import modified_within_days

from bcbio.qc import FlowcellQCMetrics, QCMetrics, SampleQCMetrics, LaneQCMetrics, FlowcellRunMetrics, SampleRunMetrics

class RunMetricsController(AbstractBaseController):
    """
    This class is an implementation of the :ref:`ICommand
    <scilifelab.pm.core.command>` interface.

    Functionality for dealing with QC data.
    """
    class Meta:
        label = 'qc'
        description = "Extension for dealing with QC data"
        arguments = [
            (['flowcell'], dict(help="Flowcell directory", nargs="?", default=None)),
            (['analysis'], dict(help="Root path to analysis folder", default=None, nargs="?")),
            (['--pre_casava'], dict(help="Toggle casava structure", default=False, action="store_true")),
            (['--project'], dict(help="Project id", default=None, action="store", type=str)),
            (['--sample'], dict(help="Sample id", default=None, action="store", type=str)),
            (['--mtime'], dict(help="Last modification time of directory (days): skip if older. Defaults to 1 day.", default=1, action="store", type=int))
            ]

    @controller.expose(hide=True)
    def default(self):
        print self._help_text


    ## Assuming standard directory structure
    def _fc_id(self):
        """Return fc id"""
        pattern = "[0-9]+_[0-9A-Za-z]+_[0-9]+_[A-Z]([A-Z0-9]+)"
        m = re.search(pattern, self.pargs.flowcell)
        return m.group(1)

    def _fc_fullname(self):
        """Return fc name (fc_date_fc_name)"""
        pattern = "([0-9]+)_[0-9A-Za-z]+_[0-9]+_([A-Z0-9]+)"
        m = re.search(pattern, self.pargs.flowcell)
        return "{}_{}".format(m.group(1), m.group(2))
    
    def _fc_parts(self):
        """return fc_name and fc_date"""
        pattern = "([0-9]+)_[0-9A-Za-z]+_[0-9]+_([A-Z0-9]+)"
        m = re.search(pattern, self.pargs.flowcell)
        return (m.group(1), m.group(2))

    ##############################
    ## Old file collection functions
    ## Will be removed when new statusdb structure in place
    ##############################
    def _collect_casava_qc_old(self):
        qc_objects = []
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(self._fc_id()))
        try:
            with open(runinfo_csv) as fh:
                runinfo_reader = csv.reader(fh)
                runinfo = [x for x in runinfo_reader]
        except IOError as e:
            self.app.log.warn(str(e))
            raise e

        for sample in runinfo[1:]:
            d = dict(zip(runinfo[0], sample))
            if self.app.pargs.project and self.app.pargs.project != d['SampleProject']:
                continue
            if self.app.pargs.sample and self.app.pargs.sample != d['SampleID']:
                continue
                
            sampledir = os.path.join(os.path.abspath(self.pargs.analysis), d['SampleProject'].replace("__", "."), d['SampleID'])
            if not os.path.exists(sampledir):
                self.app.log.warn("No such sample directory: {}".format(sampledir))
                raise IOError(2, "No such sample directory: {}".format(sampledir), sampledir)
            sample_fcdir = os.path.join(sampledir, self._fc_fullname())
            (fc_date, fc_name) = self._fc_parts()
            runinfo_yaml_file = os.path.join(sample_fcdir, "{}-bcbb-config.yaml".format(d['SampleID']))
            if not os.path.exists(runinfo_yaml_file):
                self.app.log.warn("No such yaml file for sample: {}".format(runinfo_yaml_file))
                raise IOError(2, "No such yaml file for sample: {}".format(sampledir), sampledir)
            with open(runinfo_yaml_file) as fh:
                runinfo_yaml = yaml.load(fh)
            sample_kw = dict(path=sample_fcdir, flowcell=fc_name, date=fc_date, lane=d['Lane'], barcode_name=d['SampleID'], sample_prj=d['SampleProject'].replace("__", "."), barcode_id=runinfo_yaml['details'][0]['multiplex'][0]['barcode_id'], sequence=runinfo_yaml['details'][0]['multiplex'][0]['sequence'])
            obj = SampleQCMetrics(**sample_kw)
            obj.read_picard_metrics()
            obj.parse_fastq_screen()
            obj.parse_bc_metrics()
            obj.read_fastqc_metrics()
            qc_objects.append(obj)
            
        fcdir = os.path.join(os.path.abspath(self.pargs.analysis), self.pargs.flowcell)
        fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
        fcobj = FlowcellQCMetrics(**fc_kw)
        fcobj.parse_illumina_metrics(fullRTA=False)
        qc_objects.append(fcobj)
        return qc_objects

    ## FIX ME: no way to do dry run
    def _save_obj(self, db, obj, url):
        dbobj = db.get(obj.get_db_id())
        if dbobj is None:
            obj["creation_time"] = datetime.now().isoformat()
            obj["modification_time"] = datetime.now().isoformat()
            self.app.log.info("Creating entity type %s with id %s in url %s" % (obj["entity_type"], obj.get_db_id(), url))
            db.save(obj, _id=obj["_id"])
        else:
            obj["_rev"] = dbobj.get("_rev")
            obj["creation_time"] = dbobj["creation_time"]
            ## FIXME: always ne: probably some date field gets updated somewhere
            if obj != dbobj:
                obj["modification_time"] = datetime.now().isoformat()
                self.app.log.info("Updating %s object with id %s in url %s" % (obj["entity_type"], obj.get_db_id(), url))
                db.save(obj)
            else:
                self.app.log.info("Object %s already present in url %s and not in need of updating" % (obj.get_db_id(), url))
        return True

    @controller.expose(help="OBSOLETE: upload QC metrics to couchdb")
    def upload_qc_old(self):
        if not self._check_pargs(['flowcell', 'analysis', 'url']):
            return
        if os.path.exists(os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(self._fc_id()))):
            self.log.info("Assuming casava based file structure for {}".format(self._fc_id()))
            qc_objects = self._collect_casava_qc_old()
        else:
            self.log.info("Assuming pre-casava based file structure for {}".format(self._fc_id()))
            #self._collect_pre_casava_qc_old()

        try:
            statusdb_url="http://{}:{}".format(self.pargs.url, self.pargs.port)
            couch = couchdb.Server(url=statusdb_url)
        except:
            self.app.log.warn("Connecting to server at {} failed" % statusdb_url)
            return
        self.app.log.info("Connecting to server at {} succeeded".format(statusdb_url))
        db=couch['qc']

        for obj in qc_objects:
            if self.app.pargs.debug:
                self.log.info(obj)
            else:
                self._save_obj(db, obj, statusdb_url)


    ##############################
    ## New structures
    ##############################
    def _collect_pre_casava_qc(self):
        qc_objects = []
        runinfo_yaml = os.path.join(os.path.abspath(self.pargs.flowcell), "run_info.yaml")
        try:
            with open(runinfo_yaml) as fh:
                runinfo = yaml.load(fh)
        except IOError as e:
            self.app.log.warn(str(e))
            raise e
        fcdir = os.path.abspath(self.pargs.flowcell)
        (fc_date, fc_name) = self._fc_parts()
        ## Check modification time
        if modified_within_days(fcdir, self.pargs.mtime):
            fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
            fcobj = FlowcellRunMetrics(**fc_kw)
            fcobj.parse_illumina_metrics(fullRTA=False)
            fcobj.parse_bc_metrics()
            fcobj.parse_filter_metrics()
            fcobj.parse_samplesheet_csv()
            fcobj.parse_run_info_yaml()
            qc_objects.append(fcobj)
        for info in runinfo:
            for sample in info["multiplex"]:
                sample.update({k: info.get(k, None) for k in ('analysis', 'description', 'flowcell_id', 'lane')})
                sample_kw = dict(path=fcdir, flowcell=fc_name, date=fc_date, lane=sample['lane'], barcode_name=sample['name'], sample_prj=sample['sample_prj'],
                                 barcode_id=sample['barcode_id'], sequence=sample['sequence'])
                obj = SampleRunMetrics(**sample_kw)
                obj.read_picard_metrics()
                obj.parse_fastq_screen()
                obj.parse_bc_metrics()
                obj.read_fastqc_metrics()
                qc_objects.append(obj)
        return qc_objects

    def _collect_casava_qc(self):
        qc_objects = []
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(self._fc_id()))
        try:
            with open(runinfo_csv) as fh:
                runinfo_reader = csv.reader(fh)
                runinfo = [x for x in runinfo_reader]
        except IOError as e:
            self.app.log.warn(str(e))
            raise e
        fcdir = os.path.join(os.path.abspath(self.pargs.analysis), self.pargs.flowcell)
        (fc_date, fc_name) = self._fc_parts()
        ## Check modification time
        if modified_within_days(fcdir, self.pargs.mtime):
            fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
            fcobj = FlowcellRunMetrics(**fc_kw)
            fcobj.parse_illumina_metrics(fullRTA=False)
            fcobj.parse_bc_metrics()
            fcobj.parse_samplesheet_csv()
            fcobj.parse_run_info_yaml()
            qc_objects.append(fcobj)

        for sample in runinfo[1:]:
            d = dict(zip(runinfo[0], sample))
            if self.app.pargs.project and self.app.pargs.project != d['SampleProject']:
                continue
            if self.app.pargs.sample and self.app.pargs.sample != d['SampleID']:
                continue
                
            sampledir = os.path.join(os.path.abspath(self.pargs.analysis), d['SampleProject'].replace("__", "."), d['SampleID'])
            if not os.path.exists(sampledir):
                self.app.log.warn("No such sample directory: {}".format(sampledir))
                raise IOError(2, "No such sample directory: {}".format(sampledir), sampledir)

            sample_fcdir = os.path.join(sampledir, self._fc_fullname())
            if not modified_within_days(sample_fcdir, self.pargs.mtime):
                continue
            runinfo_yaml_file = os.path.join(sample_fcdir, "{}-bcbb-config.yaml".format(d['SampleID']))
            if not os.path.exists(runinfo_yaml_file):
                self.app.log.warn("No such yaml file for sample: {}".format(runinfo_yaml_file))
                raise IOError(2, "No such yaml file for sample: {}".format(runinfo_yaml_file), runinfo_yaml_file)
            with open(runinfo_yaml_file) as fh:
                runinfo_yaml = yaml.load(fh)
            sample_kw = dict(path=sample_fcdir, flowcell=fc_name, date=fc_date, lane=d['Lane'], barcode_name=d['SampleID'], sample_prj=d['SampleProject'].replace("__", "."), barcode_id=runinfo_yaml['details'][0]['multiplex'][0]['barcode_id'], sequence=runinfo_yaml['details'][0]['multiplex'][0]['sequence'])
            obj = SampleRunMetrics(**sample_kw)
            obj.read_picard_metrics()
            obj.parse_fastq_screen()
            obj.parse_bc_metrics()
            obj.read_fastqc_metrics()
            qc_objects.append(obj)
        return qc_objects

    @controller.expose(help="Upload run metrics to statusdb")
    def upload_qc(self):
        if not self._check_pargs(['flowcell', 'analysis', 'url']):
            return
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(self._fc_id()))
        runinfo_yaml = os.path.join(os.path.abspath(self.pargs.flowcell), "run_info.yaml")
        if os.path.exists(runinfo_yaml):
            self.log.info("Assuming pre-casava based file structure for {}".format(self._fc_id()))
            qc_objects = self._collect_pre_casava_qc()
        else:
            self.log.info("Assuming casava based file structure for {}".format(self._fc_id()))
            qc_objects = self._collect_casava_qc()
        if len(qc_objects) == 0:
            self.log.info("No out-of-date qc objects for {}".format(self._fc_id()))
            return
        else:
            self.log.info("Retrieved {} updated qc objects".format(len(qc_objects)))

        ## Make sure couchdb handler is set
        if not '--couchdb' in self.app._meta.argv:
            self.app._meta.cmd_handler = 'couchdb'
            self.app._setup_cmd_handler()
        self.app.cmd.connect(self.pargs.url, self.pargs.port)
        db = self.app.cmd.db("qc")
        if not db and not self.pargs.dry_run:
            return
        for obj in qc_objects:
            if self.app.pargs.debug:
                self.log.info(obj)
                continue
            if isinstance(obj, FlowcellRunMetrics):
                self.app.cmd.save("flowcells", obj, update_fn)
            if isinstance(obj, SampleRunMetrics):
                self.app.cmd.save("samples", obj, update_fn)

def update_fn(db, obj):
    t_utc = utc_time()
    def equal(a, b):
        a_keys = [str(x) for x in a.keys() if x not in ["_id", "_rev", "creation_time", "modification_time"]]
        b_keys = [str(x) for x in b.keys() if x not in ["_id", "_rev", "creation_time", "modification_time"]]
        keys = list(set(a_keys + b_keys))
        return {k:a.get(k, None) for k in keys} == {k:b.get(k, None) for k in keys}

    if isinstance(obj, FlowcellRunMetrics):
        view = db.view("names/id_to_name")
    if isinstance(obj, SampleRunMetrics):
        view = db.view("names/id_to_name")

    d_view = {k.value:k for k in view}
    dbid =  d_view.get(obj["name"], None)
    dbobj = None
    if dbid:
        dbobj = db.get(dbid.id, None)
    if dbobj is None:
        obj["creation_time"] = t_utc
        return obj
    if equal(obj, dbobj):
        return None
    else:
        obj["creation_time"] = dbobj.get("creation_time")
        obj["modification_time"] = t_utc
        obj["_rev"] = dbobj.get("_rev")
        obj["_id"] = dbobj.get("_id")
        return obj

def load():
    """Called by the framework when the extension is 'loaded'."""
    handler.register(RunMetricsController)
