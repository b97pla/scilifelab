"""
statusdb.py
"""

from scilifelab.db import Couch

class SampleRunMetricsConnection(Couch):
    def __init__(self, **kwargs):
        super(SampleRunMetricsConnection, self).__init__(**kwargs)
        self.db = self.con["samples"]
        self.name_view = self.db.view("names/name")

    def get_flowcell(sample_name):
        """Retrieve flowcell from db for a given sample name"""
        if self.db.get(sample_name, None):
            return self.db[sample_name].get("flowcell", None)
        return None

    def get_samples(fc_id=None, sample_prj=None):
        """Retrieve samples for a given flowcell id and sample project.
        
        FIX ME: view needed with reduce operation (self.samples_view above)
        """
        pass

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

