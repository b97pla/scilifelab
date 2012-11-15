"""Database backend for connecting to statusdb"""
import re
from itertools import izip
from scilifelab.db import Couch
from scilifelab.utils.timestamp import utc_time
from uuid import uuid4

def calc_avg_qv(srm):
    """Calculate average quality score for a sample based on
    FastQC results.
    
    FastQC reports QV results in the field 'Per sequence quality scores', 
    where the subfields 'Count' and 'Quality' refer to the counts of a given quality value.
    
    :param srm: sample run metrics
    
    :returns avg_qv: Average quality value score.
    """
    try:
        count = [float(x) for x in srm["fastqc"]["stats"]["Per sequence quality scores"]["Count"]]
        quality = srm["fastqc"]["stats"]["Per sequence quality scores"]["Quality"]
        return round(sum([x*int(y) for x,y in izip(count, quality)])/sum(count), 1)
    except:
        return None

def _match_project_name_to_barcode_name(project_sample_name, sample_run_name):
    """Name mapping from project summary sample id to run info sample id"""
    if not project_sample_name.startswith("P"):
        sample_id = re.search("(\d+)_?([A-Z])?_",sample_run_name)
        if str(sample_id.group(1)) == str(project_sample_name):
            return True
        m = re.search("(_index[0-9]+)", sample_run_name)
        if not m:
            index = ""
        else:
            index = m.group(1)
        sample_id = re.search("([A-Za-z0-9\_]+)(\_index[0-9]+)?", sample_run_name.replace(index, ""))
        if str(sample_id.group(1)) == str(project_sample_name):
            return True
        else:
            return False
    if str(sample_run_name).startswith(str(project_sample_name)):
        return True
    elif str(sample_run_name).startswith(str(project_sample_name).rstrip("F")):
        return True
    elif str(sample_run_name).startswith(str(project_sample_name).rstrip("B")):
        return True
    # Add cases here
    return False

def sample_map_fn_id(sample_run_name, prj_sample):
    if 'sample_run_metrics' in prj_sample.keys():
        return prj_sample.get('sample_run_metrics').get(sample_run_name, None)
    else:
        return None

def _prune_ps_map(ps_map):
    """Only use srm_ids that end with [ACGT]+ or NoIndex
    """
    if not ps_map:
        return None
    ret = {}
    for k, v in ps_map.items():
        if re.search("_[ACGT]+$|_NoIndex$", k):
            ret[k] = v
    return ret

class SampleRunMetricsConnection(Couch):
    # FIXME: set time limits on which entries to include?
    def __init__(self, **kwargs):
        super(SampleRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["samples"]
        self.name_view = {k.key:k.id for k in self.db.view("names/name", reduce=False)}
        self.name_fc_view = {k.key:k for k in self.db.view("names/name_fc", reduce=False)}
        self.name_proj_view = {k.key:k for k in self.db.view("names/name_proj", reduce=False)}
        self.name_fc_proj_view = {k.key:k for k in self.db.view("names/name_fc_proj", reduce=False)}
    
    def _setup_views(self):
        """ """
        pass

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.

        :param name: unique name
        :param field: database field

        :returns: value if entry exists, None otherwise
        """
        self.log.debug("retrieving entry in field '{}' for name '{}'".format(field, name))
        if self.name_view.get(name, None) is None:
            self.log.warn("no field '{}' for name '{}'".format(field, name))
            return None
        if field:
            return self.db.get(self.name_view.get(name))[field]
        else:
            return self.db.get(self.name_view.get(name))

    def get_sample_ids(self, fc_id=None, sample_prj=None):
        """Retrieve sample ids subset by fc_id and/or sample_prj

        :param fc_id: flowcell id
        :param sample_prj: sample project name

        :returns sample_ids: list of couchdb sample ids
        """
        self.log.debug("retrieving sample ids subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
        fc_sample_ids = [self.name_fc_view[k].id for k in self.name_fc_view.keys() if self.name_fc_view[k].value == fc_id] if fc_id else []
        prj_sample_ids = [self.name_proj_view[k].id for k in self.name_proj_view.keys() if self.name_proj_view[k].value == sample_prj] if sample_prj else []
        # | -> union, & -> intersection
        if len(fc_sample_ids) > 0 and len(prj_sample_ids) > 0:
            sample_ids = list(set(fc_sample_ids) & set(prj_sample_ids))
        else:
            sample_ids = list(set(fc_sample_ids) | set(prj_sample_ids))
        self.log.debug("Number of samples: {}, number of fc samples: {}, number of project samples: {}".format(len(sample_ids), len(fc_sample_ids), len(prj_sample_ids)))
        return sample_ids

    def get_samples(self, fc_id=None, sample_prj=None):
        """Retrieve samples subset by fc_id and/or sample_prj

        :param fc_id: flowcell id
        :param sample_prj: sample project name

        :returns samples: list of samples
        """
        self.log.debug("retrieving samples subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
        sample_ids = self.get_sample_ids(fc_id, sample_prj)
        return [self.db.get(x) for x in sample_ids]
        
    def set_db(self):
        """Make sure we don't change db from samples"""
        pass


class FlowcellRunMetricsConnection(Couch):
    def __init__(self, **kwargs):
        super(FlowcellRunMetricsConnection, self).__init__(**kwargs)
        if not self.con:
            return
        self.db = self.con["flowcells"]
        self.name_view = {k.key:k.id for k in self.db.view("names/name", reduce=False)}

    def set_db(self):
        """Make sure we don't change db from flowcells"""
        pass

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.
        """
        self.log.debug("retrieving field entry in field '{}' for name '{}'".format(field, name))
        if self.name_view.get(name, None) is None:
            self.log.warn("no field '{}' for name '{}'".format(field, name))
            return None
        if field:
            return self.db.get(self.name_view.get(name))[field]
        else:
            return self.db.get(self.name_view.get(name))

    def get_phix_error_rate(self, name, lane, avg=True):
        """Get phix error rate"""
        fc = self.get_entry(name)
        phix_r1 = float(fc['illumina']['Summary']['read1'][lane]['ErrRatePhiX']) 
        phix_r2 = float(fc['illumina']['Summary']['read3'][lane]['ErrRatePhiX'])
        if avg:
            return (phix_r1 + phix_r2)/2
        else:
            return (phix_r1, phix_r2)/2

class ProjectSummaryConnection(Couch):
    def __init__(self, **kwargs):
        super(ProjectSummaryConnection, self).__init__(**kwargs)
        if not self.con:
            return
        self.db = self.con["projects"]
        self.name_view = {k.key:k.id for k in self.db.view("project/project_id", reduce=False)}

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.

        :param name: unique name
        :param field: database field

        :returns: value if entry exists, None otherwise
        """
        self.log.debug("retrieving field entry in field '{}' for name '{}'".format(field, name))
        if self.name_view.get(name, None) is None:
            self.log.warn("no field '{}' for name '{}'".format(field, name))
            return None
        if field:
            return self.db.get(self.name_view.get(name))[field]
        else:
            return self.db.get(self.name_view.get(name))

    def set_db(self):
        """Make sure we don't change db from projects"""
        pass

    def get_project_sample(self, project_id, sample_run_name, use_ps_map=True, use_bc_map=False,  check_consistency=False):
        """Get project sample name for a SampleRunMetrics barcode_name.
        
        :param project_id: the project id
        :param barcode_name: the barcode name of a sample run
        :param use_ps_map: use and give precedence to the project summary mapping to sample run metrics (default True)
        :param use_bc_map: use and give precedence to the barcode match mapping to sample run metrics (default False)
        :param check_consistency: use both mappings and check consistency (default False)

        :returns: project sample or None
        """
        # library_prep = re.search("P[0-9]+_[0-9]+([A-F])", sample_run_name)
        project = self.get_entry(project_id)
        if not project:
            return None
        project_samples = project.get('samples', None)
        if sample_run_name in project_samples.keys():
            return project_samples[sample_run_name]
        for project_sample_name in project_samples.keys():
            if not re.search("^P([0-9][0-9][0-9])", sample_run_name):
                sample_id = re.search("(\d+)_?([A-Z])?_",sample_run_name)
                if str(sample_id.group(1)) == str(project_sample_name):
                    return project_samples[project_sample_name]
                m = re.search("(_index[0-9]+)", sample_run_name)
                if not m:
                    index = ""
                else:
                    index = m.group(1)
                    sample_id = re.search("([A-Za-z0-9\_]+)(\_index[0-9]+)?", sample_run_name.replace(index, ""))
                    if str(sample_id.group(1)) == str(project_sample_name):
                        return project_samples[project_sample_name]
                    if str(sample_id.group(1)) == str(project_samples[project_sample_name].get("customer_name", None)):
                        return project_samples[project_sample_name]
                    # customer well names contain a 0, as in 11A07; run names don't always
                    # FIXME: a function should convert customer name to standard forms in cases like these
                    if str(sample_id.group(1)) == str(project_samples[project_sample_name].get("customer_name", None).replace("0", "")):
                        return project_samples[project_sample_name]
            else:
                if str(sample_run_name).startswith(str(project_sample_name)):
                    return project_samples[project_sample_name]
                elif str(sample_run_name).startswith(str(project_sample_name).rstrip("F")):
                    return project_samples[project_sample_name]
                elif str(sample_run_name).startswith(str(project_sample_name).rstrip("B")):
                    return project_samples[project_sample_name]
                elif str(sample_run_name).startswith(str(project_sample_name).rstrip("C")):
                    return project_samples[project_sample_name]
                elif str(sample_run_name).startswith(str(project_sample_name).rstrip("D")):
                    return project_samples[project_sample_name]
                elif str(sample_run_name).startswith(str(project_sample_name).rstrip("E")):
                    return project_samples[project_sample_name]
        return None

    def map_srm_to_name(self, project_id, include_all=True, **args):
        """Map sample run metrics names to project sample names for a
        project, possibly subset by flowcell id.

        :param project_id: project id
         :param **kw: keyword arguments to be passed to map_name_to_srm
        """
        samples = self.map_name_to_srm(project_id, **args)
        srm_to_name = {}
        for k, v in samples.items():
            if not v:
                if not include_all:
                    continue
                srm_to_name.update({"NOSRM_{}".format(k):{"sample":k, "id":None}})
            else:
                srm_to_name.update({x:{"sample":k,"id":y} for x,y in v.items()})
        return srm_to_name

    def map_name_to_srm(self, project_id, fc_id=None, use_ps_map=True, use_bc_map=False,  check_consistency=False):
        """Map project sample names to sample run metrics names for a
        project, possibly subset by flowcell id.

        :param project_id: project id
        :param fc_id: flowcell id
        :param use_ps_map: use and give precedence to the project summary mapping to sample run metrics (default True)
        :param use_bc_map: use and give precedence to the barcode match mapping to sample run metrics (default False)
        :param check_consistency: use both mappings and check consistency (default False)
        """
        project = self.get_entry(project_id)
        if not project:
            return None
        project_samples = project.get('samples', None)
        if project_samples is None:
            return None
        sample_map = {}
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        srm_samples = s_con.get_samples(fc_id=fc_id, sample_prj=project_id)
        for k, v in project_samples.items():
            sample_map[k] = None
            if check_consistency:
                use_ps_map = True
                use_bc_map = True
            if use_ps_map:
                ps_map = self._get_sample_run_metrics(v)
                sample_map[k] = _prune_ps_map(ps_map)
            if use_bc_map or not sample_map[k]:
                if not sample_map[k]: self.log.info("Using barcode map since no information in project summary for sample '{}'".format(k))
                bc_map = {s["name"]:s["_id"] for s in srm_samples if _match_project_name_to_barcode_name(k, s.get("barcode_name", None))}
                sample_map[k] = bc_map
            if check_consistency:
                if ps_map == bc_map:
                    self.log.debug("Sample {} has consistent mappings!".format(k))
                else:
                    self.log.warn("Sample {} has inconsistent mappings: ps_map {} vs barcode_match {}".format(k, ps_map, bc_map))
        return sample_map

    def _get_sample_run_metrics(self, v):
        if v.get('library_prep', None):
            library_preps = v.get('library_prep')
            return {k:v for kk in library_preps.keys() for k, v in library_preps[kk]['sample_run_metrics'].items()} if library_preps else None
        else:
            return v.get('sample_run_metrics', None)

    def get_ordered_amount(self, project_id, rounded=True, dec=1):
        """Get (rounded) ordered amount of reads in millions. 

        :param project_id: project id
        :param rounded: <boolean>
        :param dec: <integer>, number of decimal places

        :returns: ordered amount of reads if present, None otherwise
        """
        amount = self.get_entry(project_id, 'min_m_reads_per_sample_ordered')
        self.log.debug("got amount {}".format(amount))
        if not amount:
            return None
        else:
            return round(amount, dec)

    def get_qc_data(self, sample_prj, fc_id=None):
        """Get qc data for a project, possibly subset by flowcell"""
        project = self.get_entry(sample_prj)
        application = project.get("application", None) if project else None
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        samples = s_con.get_samples(fc_id=fc_id, sample_prj=sample_prj)
        qcdata = {}
        for s in samples:
            qcdata[s["name"]]={"sample":s.get("barcode_name", None),
                               "project":s.get("sample_prj", None),
                               "lane":s.get("lane", None),
                               "flowcell":s.get("flowcell", None),
                               "date":s.get("date", None),
                               "application":application,
                               "TOTAL_READS":int(s.get("picard_metrics", {}).get("AL_PAIR", {}).get("TOTAL_READS", -1)),
                               "PERCENT_DUPLICATION":s.get("picard_metrics", {}).get("DUP_metrics", {}).get("PERCENT_DUPLICATION", "-1.0"),
                               "MEAN_INSERT_SIZE":float(s.get("picard_metrics", {}).get("INS_metrics", {}).get("MEAN_INSERT_SIZE", "-1.0").replace(",", ".")),
                               "GENOME_SIZE":int(s.get("picard_metrics", {}).get("HS_metrics", {}).get("GENOME_SIZE", -1)),
                               "FOLD_ENRICHMENT":float(s.get("picard_metrics", {}).get("HS_metrics", {}).get("FOLD_ENRICHMENT", "-1.0").replace(",", ".")),
                               "PCT_USABLE_BASES_ON_TARGET":s.get("picard_metrics", {}).get("HS_metrics", {}).get("PCT_USABLE_BASES_ON_TARGET", "-1.0"),
                               "PCT_TARGET_BASES_10X":s.get("picard_metrics", {}).get("HS_metrics", {}).get("PCT_TARGET_BASES_10X", "-1.0"),
                               "PCT_PF_READS_ALIGNED":s.get("picard_metrics", {}).get("AL_PAIR", {}).get("PCT_PF_READS_ALIGNED", "-1.0"),
                               }
            target_territory = float(s.get("picard_metrics", {}).get("HS_metrics", {}).get("TARGET_TERRITORY", -1))
            pct_labels = ["PERCENT_DUPLICATION", "PCT_USABLE_BASES_ON_TARGET", "PCT_TARGET_BASES_10X",
                          "PCT_PF_READS_ALIGNED"]
            for l in pct_labels:
                if qcdata[s["name"]][l]:
                    qcdata[s["name"]][l] = float(qcdata[s["name"]][l].replace(",", ".")) * 100
            if qcdata[s["name"]]["FOLD_ENRICHMENT"] and qcdata[s["name"]]["GENOME_SIZE"] and target_territory:
                qcdata[s["name"]]["PERCENT_ON_TARGET"] = float(qcdata[s["name"]]["FOLD_ENRICHMENT"]/ (float(qcdata[s["name"]]["GENOME_SIZE"]) / float(target_territory))) * 100
        return qcdata

##############################
# Documents
##############################
class StatusDocument(dict):
    """Generic StatusDocument class"""
    _entity_type = "status_document"
    def __init__(self, **kw):
        self["_id"] = kw.get("_id", uuid4().hex)
        self["entity_type"] = self._entity_type
        self["name"] = kw.get("name", None)
        self["creation_time"] = kw.get("creation_time", utc_time())
        self["modification_time"] = None

    def __repr__(self):
        return "<{} {}>".format(self["entity_type"], self["name"])

    
class FlowcellRunMetrics(StatusDocument):
    """Flowcell level class for holding qc data."""
    _entity_type = "flowcell_run_metrics"
    def __init__(self, fc_name, fc_date, runinfo="RunInfo.xml", **kw):#, parse=True, fullRTA=False):
        StatusDocument.__init__(self, **kw)
        self["name"] = kw.get("name", "{}_{}".format(fc_name, fc_date))
        self["RunInfo"] = {"Id" : self["name"], "Flowcell":fc_name, "Date": fc_date, "Instrument": "NA"}
        self["run_info_yaml"] = kw.get("run_info_yaml", {})
        self["samplesheet_csv"] = kw.get("samplesheet_csv", [])
        self._lanes = [1,2,3,4,5,6,7,8]
        self["lanes"] = kw.get("lanes", {str(k):{"lane":str(k), "filter_metrics":{}, "bc_metrics":{}} for k in self._lanes})
        self["illumina"] = kw.get("illumina", {})

class SampleRunMetrics(StatusDocument):
    """Sample-level class for holding run metrics data"""
    _entity_type = "sample_run_metrics"
    def __init__(self, **kw):
        StatusDocument.__init__(self)
        self["entity_type"] = self._entity_type
        self["barcode_id"] = kw.get("barcode_id", None)
        self["barcode_name"] = kw.get("barcode_name", None)
        self["barcode_type"] = kw.get("barcode_type", None)
        self["bc_count"] = kw.get("bc_count", None)
        self["date"] = kw.get("date", None)
        self["flowcell"] = kw.get("flowcell", None)
        self["lane"] = kw.get("lane", None)
        self["sample_prj"] = kw.get("sample_prj", None)
        self["sequence"] = kw.get("sequence", None)
        self["barcode_type"] = kw.get("barcode_type", None)
        self["genomes_filter_out"] = kw.get("genomes_filter_out", None)
        self["name"] = "{}_{}_{}_{}".format(self["lane"], self["date"], self["flowcell"], self["sequence"])
        self["project_sample_name"] = kw.get("project_sample_name", None)
        self["project_id"] = kw.get("project_id", None)
        
        ## Metrics
        self["fastqc"] = kw.get("fastqc", {})
        self["fastq_scr"] = kw.get("fastq_scr", {})
        self["picard_metrics"] = kw.get("picard_metrics", {})

