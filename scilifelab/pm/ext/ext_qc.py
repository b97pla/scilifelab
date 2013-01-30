"""QC extension"""
import os
import csv
import yaml
import ast
import itertools
from collections import defaultdict

from cement.core import backend, controller, handler, hook
from scilifelab.utils.misc import query_yes_no
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.timestamp import modified_within_days
from scilifelab.bcbio.qc import FlowcellRunMetricsParser, SampleRunMetricsParser
from scilifelab.pm.bcbio.utils import validate_fc_directory_format, fc_id, fc_parts, fc_fullname
from scilifelab.db.statusdb import SampleRunMetricsConnection, FlowcellRunMetricsConnection, ProjectSummaryConnection, SampleRunMetricsDocument, FlowcellRunMetricsDocument
from scilifelab.utils.dry import dry
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

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
            (['--project_name'], dict(help="Project name, as in 'J.Doe_00_01'", default=None, action="store", type=str)),
            (['--project_id'], dict(help="Project id, as in 'P001'", default=None, action="store", type=str)),
            (['--names'], dict(help="Sample name mapping from barcode name to project name as a JSON string, as in \"{'sample_run_name':'project_run_name'}\". Mapping can also be given in a file", default=None, action="store", type=str)),
            (['--extensive_matching'], dict(help="Perform extensive barcode to project sample name matcing", default=False, action="store_true")),
            (['--project_alias'], dict(help="True project name as defined in project summary, as in 'J.Doe_00_01'.", default=None, action="store", type=str)),
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
        else:
            self.app.log.info("Trying to use extensive matching...")
            p_con = ProjectSummaryConnection(dbname=self.app.config.get("db", "projects"), **vars(self.app.pargs))
            project_name = self.pargs.sample_prj
            if self.pargs.project_alias:
                project_name = self.pargs.project_alias
            for s in samples:
                project_sample = p_con.get_project_sample(project_name, s["barcode_name"], extensive_matching=True)
                if project_sample:
                    self.app.log.info("using mapping '{} : {}'...".format(s["barcode_name"], project_sample["sample_name"]))
                    s["project_sample_name"] = project_sample["sample_name"]
                    s_con.save(s)
                
    ##############################
    ## New structures
    ##############################
    def _parse_samplesheet(self, runinfo, qc_objects, fc_date, fc_name, fcdir, as_yaml=False, demultiplex_stats=None):
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
                
                    parser = SampleRunMetricsParser(fcdir)
                    obj = SampleRunMetricsDocument(**sample_kw)
                    obj["picard_metrics"] = parser.read_picard_metrics(**sample_kw)
                    obj["fastq_scr"] = parser.parse_fastq_screen(**sample_kw)
                    obj["bc_count"] = parser.get_bc_count(**sample_kw)
                    obj["fastqc"] = parser.read_fastqc_metrics(**sample_kw)
                    qc_objects.append(obj)
        else:
            for sample in runinfo[1:]:
                d = dict(zip(runinfo[0], sample))
                if self.app.pargs.project_name and self.app.pargs.project_name != d['SampleProject']:
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
                obj = SampleRunMetricsDocument(**sample_kw)
                obj["picard_metrics"] = parser.read_picard_metrics(**sample_kw)
                obj["fastq_scr"] = parser.parse_fastq_screen(**sample_kw)
                obj["bc_count"] = parser.get_bc_count(demultiplex_stats=demultiplex_stats, **sample_kw)
                obj["fastqc"] = parser.read_fastqc_metrics(**sample_kw)
                qc_objects.append(obj)
        return qc_objects

    def _collect_pre_casava_qc(self):
        qc_objects = []
        as_yaml = False
        runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
        if not os.path.exists(runinfo_csv):
            LOG.warn("No such file {}: trying fallback SampleSheet.csv".format(runinfo_csv))
            runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "SampleSheet.csv")
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
            fcobj = FlowcellRunMetricsDocument(**fc_kw)
            fcobj["RunInfo"] = parser.parseRunInfo(**fc_kw)
            fcobj["RunParameters"] = parser.parseRunParameters(**fc_kw)
            fcobj["illumina"] = parser.parse_illumina_metrics(fullRTA=False, **fc_kw)
            fcobj["bc_metrics"] = parser.parse_bc_metrics(**fc_kw)
            fcobj["filter_metrics"] = parser.parse_filter_metrics(**fc_kw)
            fcobj["samplesheet_csv"] = parser.parse_samplesheet_csv(runinfo_csv=runinfo_csv, **fc_kw)
            fcobj["run_info_yaml"] = parser.parse_run_info_yaml(**fc_kw)
            qc_objects.append(fcobj)
        else:
            return qc_objects
        qc_objects = self._parse_samplesheet(runinfo, qc_objects, fc_date, fc_name, fcdir, as_yaml=as_yaml)
        return qc_objects

    def _collect_casava_qc(self):
        qc_objects = []
        runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "{}.csv".format(fc_id(self.pargs.flowcell)))
        if not os.path.exists(runinfo_csv):
            LOG.warn("No such file {}: trying fallback SampleSheet.csv".format(runinfo_csv))
            runinfo_csv = os.path.join(os.path.join(self._meta.root_path, self.pargs.flowcell), "SampleSheet.csv")
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
            fcobj = FlowcellRunMetricsDocument(fc_date, fc_name)
            fcobj["RunInfo"] = parser.parseRunInfo(**fc_kw)
            fcobj["RunParameters"] = parser.parseRunParameters(**fc_kw)
            fcobj["illumina"] = parser.parse_illumina_metrics(fullRTA=False, **fc_kw)
            fcobj["bc_metrics"] = parser.parse_bc_metrics(**fc_kw)
            fcobj["undemultiplexed_barcodes"] = parser.parse_undemultiplexed_barcode_metrics(**fc_kw)
            fcobj["illumina"].update({"Demultiplex_Stats" : parser.parse_demultiplex_stats_htm(**fc_kw)})
            fcobj["samplesheet_csv"] = parser.parse_samplesheet_csv(runinfo_csv=runinfo_csv, **fc_kw)
            qc_objects.append(fcobj)
        qc_objects = self._parse_samplesheet(runinfo, qc_objects, fc_date, fc_name, fcdir, demultiplex_stats=fcobj["illumina"]["Demultiplex_Stats"])
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
            if isinstance(obj, FlowcellRunMetricsDocument):
                dry("Saving object {}".format(repr(obj)), fc_con.save(obj))
            if isinstance(obj, SampleRunMetricsDocument):
                project_sample = p_con.get_project_sample(obj.get("sample_prj", None), obj.get("barcode_name", None), self.pargs.extensive_matching)
                if project_sample:
                    obj["project_sample_name"] = project_sample['sample_name']
                dry("Saving object {}".format(repr(obj)), s_con.save(obj))

    @controller.expose(help="Perform a multiplex QC")
    def multiplex_qc(self):
        
        MAX_UNDEMULTIPLEXED_INDEX_COUNT = 1000000
        EXPECTED_LANE_YIELD = 143000000
        MAX_PHIX_ERROR_RATE = 2.0
        MIN_PHIX_ERROR_RATE = 0.0
        read_pairs = True
        
        out_data = []
        
        if not self._check_pargs(['flowcell']):
            return
        url = self.pargs.url if self.pargs.url else self.app.config.get("db", "url")
        if not url:
            self.app.log.warn("Please provide a valid url: got {}".format(url))
            return
        
        # Construct the short form of the fcid
        sp = os.path.basename(self.pargs.flowcell).split("_")
        fcid = "_".join([sp[0],sp[-1]])
        
        # Get a connection to the flowcell database and fetch the corresponding document
        self.log.debug("Connecting to flowcell database".format(fcid))
        fc_con = FlowcellRunMetricsConnection(dbname=self.app.config.get("db", "flowcells"), **vars(self.app.pargs))
        self.log.debug("Fetching run metrics entry for flowcell {}".format(fcid))
        fc_doc = fc_con.get_entry(fcid)
        if not fc_doc:
            self.log.warn("Could not fetch run metrics entry for flowcell {}".format(fcid))
            return
   
        # Get the yield per sample from the Demultiplex_Stats
        self.log.debug("Getting yield for flowcell {}".format(fcid))
        sample_yield = self._get_yield_per_sample(fc_doc, read_pairs)
        
        # Get the yield per lane from the Demultiplex_Stats
        self.log.debug("Getting lane yield for flowcell {}".format(fcid))
        lane_yield = self._get_yield_per_lane(fc_doc, read_pairs)
        lanes = lane_yield.keys()
        
        # Get the number of samples in the pools from the Demultiplex_Stats
        self.log.debug("Getting lane pool sizes for flowcell {}".format(fcid))
        pool_size = self._get_pool_size(fc_doc)
        
        # Get the sample information from the csv samplesheet
        self.log.debug("Getting csv samplesheet data for flowcell {}".format(fcid))
        ssheet_samples = self._get_samplesheet_sample_data(fc_doc)
        if len(ssheet_samples) == 0: 
            self.log.warn("No samplesheet data available for flowcell {}".format(fcid))
        
        # Verify that all samples in samplesheet have reported metrics
        for id in ssheet_samples.keys():
            for key in ssheet_samples[id].keys():
                lane, index = key.split("_")
                project = ssheet_samples[id][key][0]
                if id not in sample_yield or \
                key not in sample_yield[id]: 
                    self.log.warn("Sample {} from project {} is in samplesheet but no yield was reported in " \
                                  "Demultiplex_Stats.htm for lane {} and index {}".format(id,
                                                                                          project,
                                                                                          lane,
                                                                                          index))
                    continue
                sample_yield[id][key].append('verified')
        
        # Check that all samples in Demultiplex_Stats have entries in Samplesheet
        for id in sample_yield.keys():
            for key in sample_yield[id].keys():
                lane, index = key.split("_")
                if "verified" not in sample_yield[id][key] and \
                index != "Undetermined":
                    self.log.warn("Sample {} from project {}, with index {} on lane {} is in Demultiplex_Stats " \
                                  "but no corresponding entry is present in SampleSheet".format(id,
                                                                                                sample_yield[id][key][1],
                                                                                                index,
                                                                                                lane))
                        
        # Check the PhiX error rate for each lane
        self.log.debug("Getting PhiX error rates for flowcell {}".format(fcid))
        for lane in lanes:
            status = "N/A"
            err_rate = fc_con.get_phix_error_rate(fcid,lane)
            if err_rate < 0:
                self.log.warn("Could not get PhiX error rate for lane {} on flowcell {}".format(lane,fcid))
            elif err_rate <= MIN_PHIX_ERROR_RATE or err_rate > MAX_PHIX_ERROR_RATE:
                status = "FAIL"
            else:
                status = "PASS"
            out_data.append([status,
                             "PhiX error rate",
                             lane,
                             err_rate,
                             "{} < PhiX e (%) <= {}".format(MIN_PHIX_ERROR_RATE,
                                                            MAX_PHIX_ERROR_RATE)])
                
        # Check that each lane received the minimum amount of reads
        for lane, reads in lane_yield.items():
            status = "FAIL"
            if reads >= EXPECTED_LANE_YIELD:
                status = "PASS"
            out_data.append([status,"Lane yield",lane,reads,"[Yield >= {}]".format(EXPECTED_LANE_YIELD)])
                
        # Check that all samples in the pool have received a minimum number of reads
        for id in sample_yield.keys():
            for key in sample_yield[id].keys():
                lane, index = key.split("_")
                if index == "Undetermined":
                    continue
                
                status = "FAIL"
                mplx_min = int(0.5*EXPECTED_LANE_YIELD/pool_size[lane])
                if sample_yield[id][key][0] >= mplx_min:
                    status = "PASS"
                out_data.append([status,"Sample yield",lane,sample_yield[id][key][1],id,sample_yield[id][key][0],"[Yield >= {}]".format(mplx_min)])
        
        # Check that the number of undetermined reads in each lane is below 10% of the total yield for the lane
        for lane, reads in lane_yield.items():
            status = "FAIL"
            key = "_".join([lane,"Undetermined"])
            undetermined = sum([counts.get(key,[0])[0] for counts in sample_yield.values()])
            cutoff = 0.1*reads
            if undetermined < cutoff:
                status = "PASS"
            out_data.append([status,"Index read",lane,undetermined,"[Undetermined < {}]".format(cutoff)])
        
        # Check that no overrepresented index sequence exists in undemultiplexed output
        self.log.debug("Fetching undemultiplexed barcode data for flowcell {}".format(fcid))
        undemux_data = self._get_undetermined_index_counts(fc_doc)
        if len(undemux_data) == 0:
            self.log.warn("No undemultiplexed barcode data available for flowcell {}".format(fcid))
        
        for lane, counts in undemux_data.items():
            mplx_min = int(min(MAX_UNDEMULTIPLEXED_INDEX_COUNT,
                               0.5*EXPECTED_LANE_YIELD/max(1,pool_size[lane])))
            status = "N/A"
            if len(counts) > 0:
                for i in range(len(counts)):
                    status = "FAIL"
                    if int(counts[i][0]) < mplx_min:
                        status = "PASS"
                    out_data.append([status,"Index",lane,counts[i][1],counts[i][2],counts[i][0],"[Undetermined index < {}]".format(mplx_min)])
            else:
                out_data.append([status,"Index",lane,"","",mplx_min,"-"])
                    
        self.app._output_data['stdout'].write("\n".join(["\t".join([str(r) for r in row]) for row in out_data]))

    def _get_undetermined_index_counts(self, fc_doc):
        """Get the top 10 undetermined index counts for each lane"""
        
        undetermined_indexes = defaultdict(list)
        undemux_data = fc_doc.get("undemultiplexed_barcodes",{})
        
        for lane, barcodes in undemux_data.items():
            for key, data in barcodes.items():
                if key != "undemultiplexed_barcodes":
                    continue
                for i in range(len(data["count"])):
                    undetermined_indexes[lane].append([data["count"][i],
                                                       data["sequence"][i],
                                                       data["index_name"][i]])
        return undetermined_indexes
        
    def _get_yield_per_sample(self, fc_doc, read_pairs=True):
        """
        Extract the yield per sample, keyed on SampleId and "Lane_Index"
        Returns a dictionary of dictionaries    
        """     
        
        # Get the yield for each sample, lane, index
        sample_yield = {}
        counts = fc_doc.get("illumina",{}).get("Demultiplex_Stats",{}).get("Barcode_lane_statistics",[])
        for sample in counts:
            id = sample['Sample ID']
            lane = sample['Lane']
            index = sample['Index']
            reads = int(sample['# Reads'].replace(',',''))
            if read_pairs:
                reads /= 2
            if id not in sample_yield:
                sample_yield[id] = {}
            sample_yield[id]["_".join([lane,index])] = [reads,sample['Project']]
            
        return sample_yield
    
    def _get_yield_per_lane(self, fc_doc, read_pairs=True):
        """
        Extract the yield per lane.
        Returns a dictionary with key-value pairs of lane and yield    
        """     
        
        # Get the yield for each lane
        lane_yield = defaultdict(int)
        counts = fc_doc.get("illumina",{}).get("Demultiplex_Stats",{}).get("Barcode_lane_statistics",[])
        for sample in counts:
            lane = sample['Lane']
            reads = int(sample['# Reads'].replace(',',''))
            if read_pairs:
                reads /= 2
            lane_yield[lane] += reads
            
        return lane_yield
    
    def _get_pool_size(self, fc_doc):
        """
        Extract the pool size for each lane.
        Returns a dictionary with key-value pairs of lane and size    
        """     
        
        pool_size = defaultdict(int)
        counts = fc_doc.get("illumina",{}).get("Demultiplex_Stats",{}).get("Barcode_lane_statistics",[])
        for sample in counts:
            lane = sample['Lane']
            index = sample['Index']
            if index != "Undetermined":
                pool_size[lane] += 1
                
        return pool_size
       
    def _get_samplesheet_sample_data(self, fc_doc):
        """
        Extract the sample data from the csv samplesheet into a dictionary of dictionaries,
        keyed with SampleId and "Lane_Index"
        """
        
        ssheet_data = fc_doc.get("samplesheet_csv",[])
        samples = {}
        for ssheet_entry in ssheet_data:
            id = ssheet_entry['SampleID']
            project = ssheet_entry['SampleProject']
            lane = ssheet_entry['Lane']
            index = ssheet_entry['Index']
            key = "_".join([lane,index])
            if id not in samples:
                samples[id] = {}
            samples[id][key] = [project]
        
        return samples

    @controller.expose(help="List the projects and corresponding applications on a flowcell")
    def list_projects(self):
        if not self._check_pargs(["flowcell"]):
            return
        
        url = self.pargs.url if self.pargs.url else self.app.config.get("db", "url")
        if not url:
            self.app.log.warn("Please provide a valid url: got {}".format(url))
            return
        if not validate_fc_directory_format(self.pargs.flowcell):
            self.app.log.warn("Path '{}' does not conform to bcbio flowcell directory format; aborting".format(self.pargs.flowcell))
            return
        
        out_data = [[self.pargs.flowcell]]
        s = self.pargs.flowcell.split("_")
        fcid = "_".join([s[0],s[-1]])
        
        self.log.debug("Establishing FlowcellRunMetricsConnection")
        fc_con = FlowcellRunMetricsConnection(dbname=self.app.config.get("db", "flowcells"), **vars(self.app.pargs))
        self.log.debug("Establishing ProjectSummaryConnection")
        p_con = ProjectSummaryConnection(dbname=self.app.config.get("db", "projects"), **vars(self.app.pargs))
    
        self.log.debug("Fetching flowcell metric document for flowcell {}".format(fcid))
        fc = fc_con.get_entry(fcid)
        if fc is None:
            self.log.warn("No flowcell metric document for flowcell {}".format(fcid))
            return
    
        self.log.debug("Fetching csv samplesheet data for flowcell {}".format(fcid))
        ssheet_data = self._get_samplesheet_sample_data(fc)
        if len(ssheet_data) == 0:
            self.log.warn("No csv samplesheet data for flowcell {}".format(fcid))
            return
    
        # Extract the project names
        projects = set([proj[0].replace("__",".") for data in ssheet_data.values() for proj in data.values()])
    
        # Extract application for each project
        for project in projects:    
            self.log.debug("Fetching project data document for project {}".format(project))
            pdoc = p_con.get_entry(project)
            if pdoc is None:
                self.log.warn("No project data document for project {}".format(project))
                continue
        
            application = pdoc.get("application","N/A")
            out_data.append([project,application])
        
        self.app._output_data['stdout'].write("\n".join(["\t".join([str(r) for r in row]) for row in out_data]))
               
def load():
    """Called by the framework when the extension is 'loaded'."""
    handler.register(RunMetricsController)
