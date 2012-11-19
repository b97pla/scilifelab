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
from scilifelab.bcbio.qc import FlowcellRunMetricsParser, SampleRunMetricsParser
from scilifelab.pm.bcbio.utils import validate_fc_directory_format, fc_id, fc_parts, fc_fullname
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectSummaryConnection, sample_run_metrics, flowcell_run_metrics
from scilifelab.utils.dry import dry

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
            (['--runqc'], dict(help="Root path to qc data folder", default=None, nargs="?")),
            (['--project'], dict(help="Project id", default=None, action="store", type=str)),
            (['--sample'], dict(help="Sample id", default=None, action="store", type=str)),
            (['--mtime'], dict(help="Last modification time of directory (days): skip if older. Defaults to 1 day.", default=1, action="store", type=int)),
            ]

    def _process_args(self):
        self._meta.root_path = self.app.pargs.runqc if self.app.pargs.runqc else self.app.config.get("runqc", "root")
        self._meta.production_root_path = self.app.config.get("runqc", "production") if self.app.config.has_option("runqc", "production") else self.app.config.get("runqc", "root")

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Add project ids and sample names to projects")
    def add(self):
        pass

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
            fc_kw = dict(fc_date = fc_date, fc_name=fc_name)
            parser = FlowcellRunMetricsParser(fcdir)
            fcobj = flowcell_run_metrics(**fc_kw)
            fcobj["illumina"] = parser.parse_illumina_metrics(fullRTA=False, **fc_kw)
            fcobj["bc_metrics"] = parser.parse_bc_metrics(**fc_kw)
            fcobj["filter_metrics"] = parser.parse_filter_metrics(**fc_kw)
            fcobj["samplesheet_csv"] = parser.parse_samplesheet_csv(**fc_kw)
            fcobj["run_info_yaml"] = parser.parse_run_info_yaml(**fc_kw)
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
                sample_kw = dict(flowcell=fc_name, date=fc_date, lane=sample['lane'], barcode_name=sample['name'], sample_prj=sample.get('sample_prj', None),
                                 barcode_id=sample['barcode_id'], sequence=sample.get('sequence', "NoIndex"))
                
                parser = SampleRunMetricsParser(fcid)
                obj = sample_run_metrics(**sample_kw)
                obj["picard_metrics"] = parser.read_picard_metrics(**sample_kw)
                obj["fastq_scr"] = parser.parse_fastq_screen(**sample_kw)
                obj["bc_metrics"] = parser.parse_bc_metrics(**sample_kw)
                obj["fastqc"] = parser.read_fastqc_metrics(**sample_kw)
                qc_objects.append(obj)
        return qc_objects

    def _collect_casava_qc(self):
        qc_objects = []
        runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
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
            fc_kw = dict(fc_date = fc_date, fc_name=fc_name)
            parser = FlowcellRunMetricsParser(fcdir)
            fcobj = flowcell_run_metrics(fc_date, fc_name)
            fcobj["illumina"] = parser.parse_illumina_metrics(fullRTA=False, **fc_kw)
            fcobj["bc_metrics"] = parser.parse_bc_metrics(**fc_kw)
            fcobj["illumina"].update({"Demultiplex_Stats" : parser.parse_demultiplex_stats_htm(**fc_kw)})
            fcobj["samplesheet_csv"] = parser.parse_samplesheet_csv(**fc_kw)
            qc_objects.append(fcobj)

        for sample in runinfo[1:]:
            d = dict(zip(runinfo[0], sample))
            if self.app.pargs.project and self.app.pargs.project != d['SampleProject']:
                continue
            if self.app.pargs.sample and self.app.pargs.sample != d['SampleID']:
                continue
                
            sampledir = os.path.join(os.path.abspath(self._meta.production_root_path), d['SampleProject'].replace("__", "."), d['SampleID'])
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
            sample_kw = dict(flowcell=fc_name, date=fc_date, lane=d['Lane'], barcode_name=d['SampleID'], sample_prj=d['SampleProject'].replace("__", "."), barcode_id=runinfo_yaml['details'][0]['multiplex'][0]['barcode_id'], sequence=runinfo_yaml['details'][0]['multiplex'][0]['sequence'])
            parser = SampleRunMetricsParser(sample_fcdir)
            obj = sample_run_metrics(**sample_kw)
            obj["picard_metrics"] = parser.read_picard_metrics(**sample_kw)
            obj["fastq_scr"] = parser.parse_fastq_screen(**sample_kw)
            obj["bc_metrics"] = parser.parse_bc_metrics(**sample_kw)
            obj["fastqc"] = parser.read_fastqc_metrics(**sample_kw)
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

        s_con = SampleRunMetricsConnection(dbname=self.app.config.get("db", "samples"), **vars(self.app.pargs))
        fc_con = FlowcellRunMetricsConnection(dbname=self.app.config.get("db", "flowcells"), **vars(self.app.pargs))
        p_con = ProjectSummaryConnection(dbname=self.app.config.get("db", "projects"), **vars(self.app.pargs))
        for obj in qc_objects:
            if self.app.pargs.debug:
                self.log.debug("{}: {}".format(str(obj), obj["_id"]))
            if isinstance(obj, flowcell_run_metrics):
                dry("Saving object {}".format(repr(obj)), fc_con.save(obj))
            if isinstance(obj, sample_run_metrics):
                project_sample = p_con.get_project_sample(obj.get("sample_prj", None), obj.get("barcode_name", None))
                if project_sample:
                    obj["project_sample_name"] = project_sample.keys()[0]
                print obj["sample_prj"], obj["project_id"]
                dry("Saving object {}".format(repr(obj)), s_con.save(obj))

def add_shared_couchdb_options(app):
    """
    Adds shared couchdb arguments to the argument object.
    
    :param app: The application object.
    
    """
    user = None
    url = None
    password = None
    if app.config.has_option("db", "user"):
        user = app.config.get("db", "user") 
    if app.config.has_option("db", "password"):
        password = app.config.get("db", "password") 
    if app.config.has_option("db", "url"):
        url = app.config.get("db", "url") 
    group = app.args.add_argument_group('couchdb', 'Options for couchdb connections')
    group.add_argument('--url', help="Database url (excluding http://). Default '{}'".format(url), default=url, nargs="?", type=str)
    group.add_argument('--port', help="Database port. Default 5984", nargs="?", default="5984", type=str)
    group.add_argument('--username', help="Database user. Default '{}'".format(user), nargs="?", default=user, type=str)
    group.add_argument('--password', help="Database password.", default=password, type=str)

def load():
    """Called by the framework when the extension is 'loaded'."""
    hook.register('post_setup', add_shared_couchdb_options)
    handler.register(RunMetricsController)
