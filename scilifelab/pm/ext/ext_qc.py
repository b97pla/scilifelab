"""QC extension"""
import os
import csv
import yaml
import ast
import itertools

from cement.core import backend, controller, handler, hook
from scilifelab.utils.misc import query_yes_no
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
            (['--sample'], dict(help="Sample id", default=None, action="store", type=str)),
            (['--mtime'], dict(help="Last modification time of directory (days): skip if older. Defaults to 1 day.", default=1, action="store", type=int)),
            (['--sample_prj'], dict(help="Sample project name, as in 'J.Doe_00_01'", default=None, action="store", type=str)),
            (['--project_id'], dict(help="Project identifier, as in 'J.Doe_00_01'", default=None, action="store", type=str)),
            (['--names'], dict(help="Sample name mapping from barcode name to project name as a JSON string, as in \"{'sample_run_name':'project_run_name'}\". Mapping can also be given in a file", default=None, action="store", type=str)),
            ]

    def _process_args(self):
        self._meta.root_path = self.app.pargs.runqc if self.app.pargs.runqc else self.app.config.get("runqc", "root")
        self._meta.production_root_path = self.app.config.get("runqc", "production") if self.app.config.has_option("runqc", "production") else self.app.config.get("runqc", "root")

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Update database objects with additional information. Currently supports updating project_id and project_sample_names in sample_run_metrics objects.")
    def update(self):
        if not self._check_pargs(["sample_prj"]):
            return
        url = self.pargs.url if self.pargs.url else self.app.config.get("db", "url")
        if not url:
            self.app.log.warn("Please provide a valid url: got {}".format(url))
            return

        s_con = SampleRunMetricsConnection(dbname=self.app.config.get("db", "samples"), **vars(self.app.pargs))
        samples = s_con.get_samples(sample_prj=self.pargs.sample_prj)

        if self.pargs.project_id:
            self.app.log.debug("Going to update 'project_id' to {} for sample runs with 'sample_prj' == {}".format(self.pargs.project_id, self.pargs.sample_prj))
            for s in samples:
                if not s.get("project_id", None) is None:
                    if not query_yes_no("'project_id':{} for sample {}; are you sure you want to overwrite?".format(s["project_id"], s["name"]), force=self.pargs.force):
                        continue
                s["project_id"] = self.pargs.project_id
                s_con.save(s)
        if self.pargs.names:
            self.app.log.debug("Going to update 'project_sample_name' for sample runs with 'sample_prj' == {}".format(self.pargs.sample_prj))
            if os.path.exists(self.pargs.names):
                with open(self.pargs.names) as fh:
                    names_d = json.load(fh)
            else:
                names_d= ast.literal_eval(self.pargs.names)
            samples_sort = sorted(samples, key=lambda s:s["barcode_name"])
            groups = {}
            for k, g in itertools.groupby(samples_sort, key=lambda x:x["barcode_name"]):
                groups[k] = list(g)
            for barcode_name in names_d:
                sample_list = groups.get(barcode_name, None)
                if not sample_list:
                    continue
                for s in sample_list:
                    if not s.get("project_sample_name", None) is None:
                        if not query_yes_no("'project_sample_name':{} for sample {}; are you sure you want to overwrite?".format(s["project_sample_name"], s["name"]), force=self.pargs.force):
                            continue
                    s["project_sample_name"] = names_d[barcode_name]
                    s_con.save(s)

                

    ##############################
    ## New structures
    ##############################
    def _parse_samplesheet(self, runinfo, qc_objects, fc_date, fc_name, as_yaml=False):
        """Parse samplesheet information and populate sample run metrics object"""
        if as_yaml:
            for info in runinfo:
                if not info.get("multiplex", None):
                    self.app.log.warn("No multiplex information for lane {}".format(info.get("lane")))
                    sample = {}
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
        else:
            for sample in runinfo[1:]:
                d = dict(zip(runinfo[0], sample))
                if self.app.pargs.project_id and self.app.pargs.project_id != d['SampleProject']:
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

    def _collect_pre_casava_qc(self):
        qc_objects = []
        as_yaml = False
        runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
        runinfo_yaml = os.path.join(os.path.abspath(self.pargs.flowcell), "run_info.yaml")
        try:
            if os.path.exists(runinfo_csv):
                with open(runinfo_csv) as fh:
                    runinfo_reader = csv.reader(fh)
                    runinfo = [x for x in runinfo_reader]
            else:
                as_yaml = True
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
        qc_objects = self._parse_samplesheet(runinfo, qc_objects, fc_date, fc_name, as_yaml=as_yaml)
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
            fcobj["undemultiplexed_barcodes"] = parser.parse_undemultiplexed_barcode_metrics(**fc_kw)
            fcobj["illumina"].update({"Demultiplex_Stats" : parser.parse_demultiplex_stats_htm(**fc_kw)})
            fcobj["samplesheet_csv"] = parser.parse_samplesheet_csv(**fc_kw)
            qc_objects.append(fcobj)
        qc_objects = self._parse_samplesheet(runinfo, qc_objects, fc_date, fc_name)
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
                dry("Saving object {}".format(repr(obj)), s_con.save(obj))

def load():
    """Called by the framework when the extension is 'loaded'."""
    handler.register(RunMetricsController)
