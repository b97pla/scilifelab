"""
statusdb.py
"""
from itertools import izip
from scilifelab.db import Couch
    
class SampleRunMetricsConnection(Couch):
    def __init__(self, **kwargs):
        super(SampleRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["samples"]
        ## FIXME: set time limits on which entries to include?
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
        self.log.info("retrieving entry in field '{}' for name '{}'".format(field, name))
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
        self.log.info("retrieving sample ids subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
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
        self.log.info("retrieving samples subset by flowcell '{}' and sample_prj '{}'".format(fc_id, sample_prj))
        sample_ids = self.get_sample_ids(fc_id, sample_prj)
        return [self.db.get(x) for x in sample_ids]

    def set_db(self):
        """Make sure we don't change db from samples"""
        pass

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
            return None



class FlowcellRunMetricsConnection(Couch):
    def __init__(self, **kwargs):
        super(FlowcellRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["flowcells"]
        self.name_view = {k.key:k.id for k in self.db.view("names/name", reduce=False)}

    def set_db(self):
        """Make sure we don't change db from flowcells"""
        pass

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.
        """
        self.log.info("retrieving field entry in field '{}' for name '{}'".format(field, name))
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

class ProjectSummary(Couch):
    def __init__(self, **kwargs):
        super(ProjectSummary, self).__init__(**kwargs)
        self.db = self.con["projects"]

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.

        :param name: unique name
        :param field: database field

        :returns: value if entry exists, None otherwise
        """
        self.log.info("retrieving field entry in field '{}' for name '{}'".format(field, name))
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
        self.log.info("retrieving field entry in field '{}' for name '{}'".format(field, name))
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
        
