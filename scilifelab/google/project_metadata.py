#!/usr/bin/env python
"""Module for extracting a project's UppNex id from a spreadsheet on Google Docs"""

import os 
from scilifelab.google.google_docs import SpreadSheet
from scilifelab.google import get_credentials

class ProjectMetaData:
    """A placeholder for metadata associated with a project.
    The data are fetched from the Google Docs spreadsheet.
    """

    def __init__(self, project_name, config):
        """Initialize the object"""

        # Map internal attribute names to the GPL column headers
        col_mapping = self.column_mapping()
        for attr in col_mapping.keys():
            setattr(self, attr, None)

        # Get the name of the spreadsheet where uppnex ids can be found
        gdocs_config = config.get("gdocs", config.get("gdocs_upload",{}))
        cred_file = gdocs_config.get("credentials_file",gdocs_config.get("gdocs_credentials"))
        ssheet_title = gdocs_config.get("projects_spreadsheet")
        wsheet_title = gdocs_config.get("projects_worksheet")

        # Get the credentials
        credentials = get_credentials(cred_file)
        assert credentials is not None, \
        "The Google Docs credentials could not be found."
        assert ssheet_title is not None and wsheet_title is not None, \
            "The names of the projects spreadsheet and worksheet on Google \
            Docs could not be found."

        # Connect to the spread- and worksheet
        ssheet = SpreadSheet(credentials, ssheet_title)
        assert ssheet is not None, \
            "Could not fetch '{}' from Google Docs.".format(ssheet_title)

        # We allow multiple, comma-separated worksheets to be searched
        for wtitle in wsheet_title.split(','):
            wsheet = ssheet.get_worksheet(wtitle.strip())
            if not wsheet:
                print("WARNING: Could not locate {} in {}".format(wsheet_title, ssheet_title))
                continue

            # Get the rows for the project
            rows = ssheet.get_cell_content(wsheet)
            header = ssheet.get_header(wsheet)
            column_indexes = {attr: ssheet.get_column_index(wsheet,col)-1 for attr, col in col_mapping.items()}
            for row in rows:
                # skip if this is not the project we're interested in
                if row[column_indexes["project_name"]] != project_name:
                    continue
                
                # Will only use the first result found to set each attribute
                for attr, index in column_indexes.items():
                    setattr(self, attr, row[index])

                # We have found the project data so stop iterating
                return

    @staticmethod
    def column_mapping():
        """Returns a dictionary with mapping between class attributes
        and the corresponding Google docs spreadsheet column name
        """

        return {
	       "type" : "Type",
           "ref_genome":"Ref genome",
           "project_id": "ID",
           "project_name": "Project name",
           "queue_date": "Queue date",
           "no_samples": "No of samples (in order)",
           "lanes_plates": "Lanes / Plates",
           "min_reads_per_sample": "minimal M read pairs/sample (passed filter)",
           "customer_reference": "Customer reference",
           "application": "Application",
           "no_finished_samples": "No of samples finished (All sequencing finished)",
           "uppnex_id": "Uppnex ID"
        }

    def _get_project_id(self):
        return self._project_id

    def _set_project_id(self, value):
        self._project_id = value

    def _get_project_name(self):
        return self._project_name

    def _set_project_name(self, value):
        self._project_name = value

    def _get_queue_date(self):
        return self._queue_date

    def _set_queue_date(self, value):
        self._queue_date = value

    def _get_no_samples(self):
        return self._no_samples

    def _set_no_samples(self, value):
        self._no_samples = value

    def _get_lanes_plates(self):
        return self._lanes_plates

    def _set_lanes_plates(self, value):
        self._lanes_plates = value

    def _get_min_reads_per_sample(self):
        return self._min_reads_per_sample

    def _set_min_reads_per_sample(self, value):
        self._min_reads_per_sample = value

    def _get_customer_reference(self):
        return self._customer_reference

    def _set_customer_reference(self, value):
        self._customer_reference = value

    def _get_application(self):
        return self._application

    def _set_application(self, value):
        self._application = value

    def _get_no_finished_samples(self):
        return self._no_finished_samples

    def _set_no_finished_samples(self, value):
        self._no_finished_samples = value

    def _get_uppnex_id(self):
        return self._uppnex_id

    def _set_uppnex_id(self, value):
        self._uppnex_id = value

    def _get_type(self):
        return self._type

    def _set_type(self, value):
        self._type = value

    def _get_ref_genome(self):
        return self._ref_genome
    
    def _set_ref_genome(self, value):
        self._ref_genome = value
 
    type = property(_get_type, _set_type)
    project_id = property(_get_project_id, _set_project_id)
    project_name = property(_get_project_name, _set_project_name)
    queue_date = property(_get_queue_date, _set_queue_date)
    no_samples = property(_get_no_samples, _set_no_samples)
    lanes_plates = property(_get_lanes_plates, _set_lanes_plates)
    min_reads_per_sample = property(_get_min_reads_per_sample, _set_min_reads_per_sample)
    customer_reference = property(_get_customer_reference, _set_customer_reference)
    application = property(_get_application, _set_application)
    no_finished_samples = property(_get_no_finished_samples, _set_no_finished_samples)
    uppnex_id = property(_get_uppnex_id, _set_uppnex_id)
    ref_genome = property(_get_ref_genome,_set_ref_genome)
    
