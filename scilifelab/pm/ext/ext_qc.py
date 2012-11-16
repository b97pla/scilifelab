"""QC extension"""
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
from scilifelab.bcbio.qc import FlowcellRunMetrics, SampleRunMetrics, QCLEVELS
from scilifelab.pm.bcbio.utils import validate_fc_directory_format, fc_id, fc_parts, fc_fullname

class RunMetricsController(AbstractBaseController):
    """
    This class is an implementation of the :ref:`ICommand
    <scilifelab.pm.core.command>` interface.

    Functionality for dealing with QC data.
    """
    class Meta:
        label = 'qc'
        description = "Extension for dealing with QC data"
        arguments = [(['flowcell'], dict(help="Flowcell directory", nargs="?", default=None)),
                     (['--runqc'],  dict(help="Root path to qc data folder", default=None, nargs="?")),
                     (['--project'], dict(help="Project id", default=None, action="store", type=str)),
                     (['--sample'], dict(help="Sample id", default=None, action="store", type=str)),
                     (['--mtime'], dict(help="Last modification time of directory (days): skip if older. Defaults to 1 day.", default=1, action="store", type=int))]


    def _process_args(self):
        self._meta.root_path = self.app.pargs.runqc if self.app.pargs.runqc else self.app.config.get("runqc", "root")

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

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
        (fc_date, fc_name) = fc_parts(self.pargs.flowcell)
        ## Check modification time
        if modified_within_days(fcdir, self.pargs.mtime):
            fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
            fcobj = FlowcellRunMetrics(**fc_kw)
            fcobj.parse_illumina_metrics(fullRTA=False)
            fcobj.parse_bc_metrics()
            fcobj.parse_filter_metrics()
            if not fcobj.parse_samplesheet_csv():
                fcobj.parse_run_info_yaml()
            qc_objects.append(fcobj)
        else:
            return qc_objects
        for info in runinfo:
            if not info.get("multiplex", None):
                self.app.log.warn("No multiplex information for lane {}".format(info.get("lane")))
                sample.update({k: info.get(k, None) for k in ('analysis', 'description', 'flowcell_id', 'lane')})
                sample_kw = dict(path=fcdir, flowcell=fc_name, date=fc_date, lane=sample.get('lane', None), barcode_name=sample.get('name', None), sample_prj=sample.get('sample_prj', None),
                                 barcode_id=sample.get('barcode_id', None), sequence=sample.get('sequence', "NoIndex"))
            for sample in info["multiplex"]:
                sample.update({k: info.get(k, None) for k in ('analysis', 'description', 'flowcell_id', 'lane')})
                sample_kw = dict(path=fcdir, flowcell=fc_name, date=fc_date, lane=sample['lane'], barcode_name=sample['name'], sample_prj=sample.get('sample_prj', None),
                                 barcode_id=sample['barcode_id'], sequence=sample.get('sequence', "NoIndex"))
                obj = SampleRunMetrics(**sample_kw)
                obj.read_picard_metrics()
                obj.parse_fastq_screen()
                obj.parse_bc_metrics()
                obj.read_fastqc_metrics()
                qc_objects.append(obj)
        return qc_objects

    def _collect_casava_qc(self):
        qc_objects = []
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
        try:
            with open(runinfo_csv) as fh:
                runinfo_reader = csv.reader(fh)
                runinfo = [x for x in runinfo_reader]
        except IOError as e:
            self.app.log.warn(str(e))
            raise e
        fcdir = os.path.join(os.path.abspath(self._meta.root_path), self.pargs.flowcell)
        (fc_date, fc_name) = fc_parts(self.pargs.flowcell)
        ## Check modification time
        if modified_within_days(fcdir, self.pargs.mtime):
            fc_kw = dict(path=fcdir, fc_date = fc_date, fc_name=fc_name)
            fcobj = FlowcellRunMetrics(**fc_kw)
            fcobj.parse_illumina_metrics(fullRTA=False)
            fcobj.parse_bc_metrics()
            fcobj.parse_demultiplex_stats_htm()
            fcobj.parse_samplesheet_csv()
            qc_objects.append(fcobj)

        for sample in runinfo[1:]:
            d = dict(zip(runinfo[0], sample))
            if self.app.pargs.project and self.app.pargs.project != d['SampleProject']:
                continue
            if self.app.pargs.sample and self.app.pargs.sample != d['SampleID']:
                continue
                
            sampledir = os.path.join(os.path.abspath(self._meta.root_path), d['SampleProject'].replace("__", "."), d['SampleID'])
            if not os.path.exists(sampledir):
                self.app.log.warn("No such sample directory: {}".format(sampledir))
                continue
            sample_fcdir = os.path.join(sampledir, fc_fullname(self.pargs.flowcell))
            if not os.path.exists(sample_fcdir):
                self.app.log.warn("No such sample flowcell directory: {}".format(sample_fcdir))
                continue
            if not modified_within_days(sample_fcdir, self.pargs.mtime):
                continue
            runinfo_yaml_file = os.path.join(sample_fcdir, "{}-bcbb-config.yaml".format(d['SampleID']))
            if not os.path.exists(runinfo_yaml_file):
                self.app.log.warn("No such yaml file for sample: {}".format(runinfo_yaml_file))
                raise IOError(2, "No such yaml file for sample: {}".format(runinfo_yaml_file), runinfo_yaml_file)
            with open(runinfo_yaml_file) as fh:
                runinfo_yaml = yaml.load(fh)
            if not runinfo_yaml['details'][0].get("multiplex", None):
                self.app.log.warn("No multiplex information for sample {}".format(d['SampleID']))
                continue
            sample_kw = dict(path=sample_fcdir, flowcell=fc_name, date=fc_date, lane=d['Lane'], barcode_name=d['SampleID'], sample_prj=d['SampleProject'].replace("__", "."), barcode_id=runinfo_yaml['details'][0]['multiplex'][0]['barcode_id'], sequence=runinfo_yaml['details'][0]['multiplex'][0]['sequence'])
            obj = SampleRunMetrics(**sample_kw)
            obj.read_picard_metrics()
            obj.parse_fastq_screen()
            obj.parse_bc_metrics()
            obj.read_fastqc_metrics()
            if d['Operator'].endswith("failed"):
                obj["qcstatus"]["instrument_run_qc"] = QCLEVELS["FAIL"]
            else:
                obj["qcstatus"]["instrument_run_qc"] = QCLEVELS["PASS"]
            qc_objects.append(obj)
        return qc_objects

    @controller.expose(help="Upload run metrics to statusdb")
    def upload_qc(self):
        if not self._check_pargs(['flowcell']):
            return
        url = self.pargs.url if self.pargs.url else self.app.config.get("db", "url")
        if not url:
            self.app.log.warn("Please provide a valid url: got {}".format(url))
            return
        if not validate_fc_directory_format(self.pargs.flowcell):
            self.app.log.warn("Path '{}' does not conform to bcbio flowcell directory format; aborting".format(self.pargs.flowcell))
            return
            
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
        runinfo_yaml = os.path.join(os.path.abspath(self.pargs.flowcell), "run_info.yaml")
        (fc_date, fc_name) = fc_parts(self.pargs.flowcell)
        if int(fc_date) < 120815:
            self.log.info("Assuming pre-casava based file structure for {}".format(fc_id(self.pargs.flowcell)))
            qc_objects = self._collect_pre_casava_qc()
        else:
            self.log.info("Assuming casava based file structure for {}".format(fc_id(self.pargs.flowcell)))
            qc_objects = self._collect_casava_qc()

        if len(qc_objects) == 0:
            self.log.info("No out-of-date qc objects for {}".format(fc_id(self.pargs.flowcell)))
            return
        else:
            self.log.info("Retrieved {} updated qc objects".format(len(qc_objects)))

        self.app._meta.cmd_handler = 'couchdb'
        self.app._setup_cmd_handler()
        self.app.cmd.connect(url, self.pargs.port)
        for obj in qc_objects:
            if self.app.pargs.debug:
                self.log.debug("{}: {}".format(str(obj), obj["_id"]))
                continue
            if isinstance(obj, FlowcellRunMetrics):
                self.app.cmd.save("flowcells", obj, update_fn)
            if isinstance(obj, SampleRunMetrics):
                self.app.cmd.save("samples", obj, update_fn)

def update_fn(db, obj):
    """Compare object with object in db if present.

    :param db: couch database
    :param obj: database object to save

    :returns: database object to save and database id if present
    """
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
        return (obj, dbid)
    if equal(obj, dbobj):
        return (None, dbid)
    else:
        obj["creation_time"] = dbobj.get("creation_time")
        obj["modification_time"] = t_utc
        obj["_rev"] = dbobj.get("_rev")
        obj["_id"] = dbobj.get("_id")
        return (obj, dbid)

def load():
    """Called by the framework when the extension is 'loaded'."""
    handler.register(RunMetricsController)
