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

from scilifelab.bcbio.qc import FlowcellRunMetrics, SampleRunMetrics

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
            ## FIXME: analysis is a confusing name
            (['analysis'], dict(help="Root path to analysis folder", default=None, nargs="?")),
            (['--pre_casava'], dict(help="Toggle casava structure", default=False, action="store_true")),
            (['--project'], dict(help="Project id", default=None, action="store", type=str)),
            (['--sample'], dict(help="Sample id", default=None, action="store", type=str)),
            (['--mtime'], dict(help="Last modification time of directory (days): skip if older. Defaults to 1 day.", default=1, action="store", type=int)),
            (['--remap'], dict(help="WARNING: temporary fix to remap objects to LANE_DATE_FC_BCSEQUENCE.", default=False, action="store_true")),
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
            #fcobj.parse_bc_metrics()
            fcobj.parse_demultiplex_stats_htm()
            fcobj.parse_samplesheet_csv()
            #fcobj.parse_run_info_yaml()
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
                continue
            sample_fcdir = os.path.join(sampledir, self._fc_fullname())
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
            qc_objects.append(obj)
        return qc_objects

    @controller.expose(help="Upload run metrics to statusdb")
    def upload_qc(self):
        if self.pargs.remap:
            self.app.log.info("using update_remap_fn")
            update_fn_ptr = update_remap_fn
        else:
            self.app.log.info("using update_fn")
            update_fn_ptr = update_fn
        if not self._check_pargs(['flowcell', 'analysis', 'url']):
            return
        runinfo_csv = os.path.join(os.path.abspath(self.pargs.flowcell), "{}.csv".format(self._fc_id()))
        runinfo_yaml = os.path.join(os.path.abspath(self.pargs.flowcell), "run_info.yaml")
        ## Have to set pre-casava flag!
        if self.pargs.pre_casava:
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
        for obj in qc_objects:
            if self.app.pargs.debug:
                self.log.debug("{}: {}".format(str(obj), obj["_id"]))
                continue
            if isinstance(obj, FlowcellRunMetrics):
                self.app.cmd.save("flowcells", obj, update_fn_ptr)
            if isinstance(obj, SampleRunMetrics):
                self.app.cmd.save("samples", obj, update_fn_ptr)

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

def update_remap_fn(db, obj):
    t_utc = utc_time()
    def equal(a, b):
        a_keys = [str(x) for x in a.keys() if x not in ["_id", "_rev", "creation_time", "modification_time"]]
        b_keys = [str(x) for x in b.keys() if x not in ["_id", "_rev", "creation_time", "modification_time"]]
        keys = list(set(a_keys + b_keys))
        return {k:a.get(k, None) for k in keys} == {k:b.get(k, None) for k in keys}

    if isinstance(obj, FlowcellRunMetrics):
        print "Working on FlowcellRunMetrics"
        view = db.view("names/id_to_name")
        obj_name = obj["name"]
    if isinstance(obj, SampleRunMetrics):
        print "Working on SampleRunMetrics"
        view = db.view("names/id_to_name")
        sample_bc_to_seq = {}
        sample_seq_to_bc = {}
        for k in db:
            dbobj = db.get(k)
            if not dbobj:
                continue
            if not "name" in dbobj:
                print dbobj
                
            sample_seq_id = "{}_{}_{}_{}".format(obj.get("lane"), obj.get("date"), obj.get("flowcell"), obj.get("sequence", "NoIndex"))
            if not sample_seq_id in sample_seq_to_bc.keys():
                sample_seq_to_bc[sample_seq_id] = [k]
            else:
                sample_seq_to_bc[sample_seq_id].append(k)
            if not dbobj["name"] in sample_bc_to_seq.keys():
                sample_bc_to_seq[dbobj["name"]] = [k]
            else:
                sample_bc_to_seq[dbobj["name"]].append(k)
                obj_name = "{}_{}_{}_{}".format(obj["lane"], obj["date"], obj["flowcell"], obj["barcode_id"])
            for k,v in sample_bc_to_seq.items():
                if len(v) > 1:
                    print "WARNING: duplicate sample_bcid entry: {}, {}".format(k, v)
            for k,v in sample_seq_to_bc.items():
                if len(v) > 1:
                    print "WARNING: duplicate sample_seq_id entry: {}, {}".format(k, v)

    d_view = {k.value:k for k in view}
    dbid = d_view.get(obj_name, None)
    ## Make sure that dbid only has single entry for sample run metrics
    if isinstance(obj, SampleRunMetrics):
        if len(sample_bc_to_seq[dbid.id]):
            self.log.warn("Duplicate entry for sample {}, barcode name {}, sample seq id {}".format(obj["barcode_name"], obj_bc_name, sample_bc_to_seq[dbid.id]))
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
