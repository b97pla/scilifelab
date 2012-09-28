"""
statusdb.py
"""

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

    def get_entry(self, sample_name, entry):
        """Retrieve entry from db for a given sample name.

        :param sample_name: sample name
        :param entry: database field

        :returns: value if entry exists, None otherwise
        """
        self.log.info("retrieving entry '{}' for sample '{}'".format(entry, sample_name))
        if not self.name_view.get(sample_name, None) is None:
            return self.db.get(self.name_view.get(sample_name))[entry]
        self.log.warn("no entry '{}' for sample '{}'".format(entry, sample_name))
        return None

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

class FlowcellRunMetricsConnection(Couch):
    def __init__(self, **kwargs):
        self.name_view = "names/name"
        super(SampleRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["flowcells"]


class ProjectSummary(Couch):
    def __init__(self, **kwargs):
        pass

