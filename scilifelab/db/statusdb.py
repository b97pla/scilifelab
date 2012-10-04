"""Database backend for connecting to statusdb"""
import re
from itertools import izip
from scilifelab.db import Couch

def sample_map_fn(sample_run_name, project_sample_name):
    """Name mapping from project summary sample id to run info sample id"""
    if not project_sample_name.startswith("P"):
        try:
            sid = re.search("(\d+)([A-Z])?_",sample_run_name)
            if str(sid.group(1)) == str(project_sample_name):
                return True
            else: 
                return False
        except:
            pass
    if str(sample_run_name).startswith(str(project_sample_name)):
        return True
    elif str(sample_run_name).startswith(str(project_sample_name).rstrip("F")):
        return True
    elif str(sample_run_name).startswith(str(project_sample_name).rstrip("B")):
        return True
    ## Add cases here
    return False

def sample_map_fn_id(sample_run_name, prj_sample):
    if 'sample_run_metrics' in prj_sample.keys():
        return prj_sample.get('sample_run_metrics').get(sample_run_name, None)
    else:
        return None

    
class SampleRunMetricsConnection(Couch):
    ## FIXME: set time limits on which entries to include?
    def __init__(self, **kwargs):
        super(SampleRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["samples"]
        self.name_view = {k.key:k.id for k in self.db.view("names/name", reduce=False)}
        self.name_fc_view = {k.key:k for k in self.db.view("names/name_fc", reduce=False)}
        self.name_proj_view = {k.key:k for k in self.db.view("names/name_proj", reduce=False)}
        self.name_fc_proj_view = {k.key:k for k in self.db.view("names/name_fc_proj", reduce=False)}

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

    def get_sample_ids(self, fc_id, sample_prj=None):
        """Retrieve sample ids subset by fc_id and possibly sample_prj

        :param fc_id: flowcell id
        :param sample_prj: sample project name

        :returns sample_ids: list of couchdb sample ids
        """
        self.log.debug("retrieving sample ids subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
        sample_ids = [self.name_fc_view[k].id for k in self.name_fc_view.keys() if self.name_fc_view[k].value == fc_id]
        if sample_prj:
            prj_sample_ids = [self.name_fc_view[k].id for k in self.name_proj_view.keys() if self.name_proj_view[k].value == sample_prj]
            sample_ids = list(set(sample_ids).intersection(set(prj_sample_ids)))
        return sample_ids

    def get_samples(self, fc_id, sample_prj=None):
        """Retrieve samples subset by fc_id and possibly sample_prj

        :param fc_id: flowcell id
        :param sample_prj: sample project name

        :returns samples: list of samples
        """
        self.log.debug("retrieving samples subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
        sample_ids = self.get_sample_ids(fc_id, sample_prj)
        return [self.db.get(x) for x in sample_ids]

    def get_project_sample_ids(self, sample_prj):
        """Retrieve sample ids subset by sample_prj

        :param sample_prj: sample project name

        :returns sample_ids: list of couchdb sample ids
        """
        self.log.debug("retrieving sample ids subset by sample_prj '{}'".format(sample_prj))
        sample_ids = [self.name_proj_view[k].id for k in self.name_proj_view.keys() if self.name_proj_view[k].value == sample_prj]
        return sample_ids

    def get_project_samples(self, sample_prj):
        """Retrieve samples subset related to a project.

        :param sample_prj: sample project name

        :returns samples: list of samples
        """
        self.log.debug("retrieving samples subset by sample_prj '{}'".format(sample_prj))
        sample_ids = self.get_project_sample_ids(sample_prj)
        return [self.db.get(x) for x in sample_ids]
        
    def set_db(self):
        """Make sure we don't change db from samples"""
        pass

    ## FIX ME: operations on sample run metrics objects should be
    ## separated from the connection. Either implement a
    ## sample_run_metrics object (subclassing ViewResults) with this
    ## function or move to utils or similar
    def calc_avg_qv(self, name):
        """Calculate average quality score for a sample based on
        FastQC results.
        
        FastQC reports QV results in the field 'Per sequence quality scores', 
        where the subfields 'Count' and 'Quality' refer to the counts of a given quality value.
        
        :param name: sample name
        
        :returns avg_qv: Average quality value score.
        """
        srm = self.get_entry(name)
        try:
            count = [float(x) for x in srm["fastqc"]["stats"]["Per sequence quality scores"]["Count"]]
            quality = srm["fastqc"]["stats"]["Per sequence quality scores"]["Quality"]
            return round(sum([x*int(y) for x,y in izip(count, quality)])/sum(count), 1)
        except:
            self.log.warn("Calculation of average quality failed for sample {}, id {}".format(srm["name"], srm["_id"]))
            return None

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

    def map_sample_run_to_project_name(self, project_id):
        """Map sample run to project name"""
        project = self.get_entry(project_id)
        project_samples = project.get('samples', None)
        d = {pid:project_samples[pid].get('sample_run_metrics', None) for pid in project_samples.keys()}
        retlist = []
        for pid in d.keys():
            if d[pid]:
                for sid in d[pid].keys():
                    dbid = d[pid][sid]
                    retlist.append((pid, sid, dbid))
            else:
                retlist.append((pid, None, None))
        return retlist
        
    def sample_map(self, project_id, fc_id=None):
        """Map project sample names to sample run names for a project,
        possibly subset by flowcell id.

        """
        project = self.get_entry(project_id)
        if not project:
            return None
        project_samples = project.get('samples', None)
        if project_samples is None:
            return None
        sample_map = {k:v.get('sample_run_metrics', None) for k, v in project_samples.items()}
        def check_duplicates(sid,d):
            if d is None:
                return
            keys = ["_".join(x.split("_")[0:3]) for x in d.keys()]
            if len(keys) != len(list(set(keys))):
                self.log.warn("Sample run metrics duplicate for sample {}: {}".format(sid,d))
        all_runs = {}
        for k,v in sample_map.items():
            check_duplicates(k, v)   
            if not v is None:
                for x in v.keys():
                    if not x in all_runs.keys():
                        all_runs[x] = [k]
                    else:
                        all_runs[x].append(k)
        for k, v in all_runs.items():
            if len(v) > 1:
                self.log.warn("Sample run metrics duplicate over samples: sample run metrics id {} maps to {}".format(k, v))
        return sample_map
        
    def map_sample_run_names(self, project_id, fc_id=None):
        """Map sample run names for a project, possibly subset by fc_id, to
        sample name defined by project.

        :param project_id: project id
        :param fc_id: flowcell id

        :returns: dict with key sample run name and value dict with keys 'sample_id' and 'project_sample'
        """
        project = self.get_entry(project_id)
        project_samples = project.get('samples', None)

        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        if fc_id:
            sample_run_samples = s_con.get_samples(fc_id, project_id)
        else:
            sample_run_samples = s_con.get_project_samples(project_id)

        ## FIX ME: This code is inefficient. Want mapping from
        ## sample_run_metrics name to project sample name (and vice
        ## versa) Information is (sometimes) present in
        ## project_summary object as database id, otherwise need match
        
        ## sample names: get sample_run_names from project_samples if present
        sample_map = {}
        for x in sample_run_samples:
            sample_map[x["name"]] = None
            for y in project_samples.keys():
                srm_id = sample_map_fn_id(x["name"], project_samples[str(y)])
                self.log.debug("project sample name {} is mapped to sample_run_metrics id {}".format(srm_id, x["_id"]))
                if sample_map_fn(x["barcode_name"], y):
                    self.log.debug("sample_sample_run_metrics name {}, barcode name {}, id {},  is mapped to project sample name {}".format(x["name"], x["barcode_name"], x["_id"], y))
                    sample_map[str(x["name"])] = {'sample_id':x["_id"], 'project_sample':y}
                if not srm_id is None:
                    self.log.debug("{} {}".format(srm_id, x["_id"]))
        for x in sample_map:
            if x is None:
                self.log.warn("No mapping from name to sample for {}".format(x))
        return sample_map

    def get_ordered_amount(self, project_id, rounded=True, dec=1):
        """Get (rounded) ordered amount of reads in millions. 

        :param project_id: project id
        :param rounded: <boolean>
        :param dec: <integer>, number of decimal places

        :returns: ordered amount of reads if present, None otherwise
        """
        amount = self.get_entry(project_id, 'min_m_reads_per_sample_ordered')
        if not amount:
            return None
        else:
            return round(amount, dec)

class ProjectQCSummaryConnection(Couch):
    """ProjectQCSummary. Connection to old QC database.

    """
    def __init__(self, **kwargs):
        super(ProjectQCSummaryConnection, self).__init__(**kwargs)
        self.db = self.con["qc"]
        self.name_view = {self.db.get(k.id)["Project_id"]:k.id for k in self.db.view("entitytypes/ProjectSummary", reduce=False)}

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.

        :param name: unique name
        :param field: database field

        :returns: value if entry exists, None otherwise
        """
        self.log.debug("retrieving field entry in field '{}' for name '{}'".format(field, name))
        if self.name_view.get(name, None) is None:
            self.log.warn("no entry with name '{}'".format(name))
            return None
        if field:
            return self.db.get(self.name_view.get(name))[field]
        else:
            return self.db.get(self.name_view.get(name))

    def set_db(self):
        """Make sure we don't change db from projects"""
        pass

    def get_customer_name(self, sample_name):
        """Get a customer name corresponding to <sample_name>.

        :param sample_name: sample name

        :returns: customer name if it exists, None otherwise
        """
        pass

    def map_sample_run_names(self, project_id, fc_id):
        """Map sample run names for a project subset by fc_id to
        sample name defined by project.

        :param project_id: flowcell id
        :param fc_id: flowcell id

        :returns: dict with key sample run name and value project sample name
        """
        project = self.get_entry(project_id)
        project_samples = project['Samples']
        s_con = SampleRunMetricsConnection(username=self.user, password=self.pw, url=self.url)
        sample_run_samples = s_con.get_samples(fc_id, project_id)

        ## FIX ME: Mapping must be much more general to account for anomalous cases
        def sample_map_fn(sample_run_name, project_sample_name):
            """Mapping from project summary sample id to run info sample id"""
            if str(sample_run_name).startswith(str(project_sample_name)):
                return True
            ## Add cases here
            return False

        sample_map = {}
        for x in sample_run_samples:
            sample_map[x["barcode_name"]] = None
            for y in project_samples.keys():
                if sample_map_fn(x["barcode_name"], y):
                    sample_map[str(x["barcode_name"])] = y
        return sample_map

    def get_samples(self):
        """Get samples defined in ProjectQCSummary document.

        :returns: list of sample names
        """
        pass
        
