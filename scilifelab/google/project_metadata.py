#!/usr/bin/env python
"""Module for extracting a project's UppNex id from a spreadsheet on Google Docs"""

import bcbio.google
import bcbio.google.spreadsheet
import os 


class ProjectMetaData:
    """A placeholder for metadata associated with a project.
    The data are fetched from the Google Docs spreadsheet.
    """

    def __init__(self, project_name, config):
        """Initialize the object"""

        col_mapping = self.column_mapping()
        name_column = col_mapping["project_name"]
        columns = []
        attributes = []
        for attr, col in col_mapping.items():
            columns.append(col)
            attributes.append(attr)
            setattr(self, attr, None)

        # Get the credentials
        encoded_credentials = bcbio.google.get_credentials(config)
        assert encoded_credentials is not None, \
        "The Google Docs credentials could not be found."

        # Get the name of the spreadsheet where uppnex ids can be found
        gdocs_config = config.get("gdocs_upload", {})
        ssheet_title = gdocs_config.get("projects_spreadsheet", None)
        wsheet_title = gdocs_config.get("projects_worksheet", None)

        assert ssheet_title is not None and wsheet_title is not None, \
            "The names of the projects spreadsheet and worksheet on Google \
            Docs could not be found."

        # Connect to the spread- and worksheet
        client = bcbio.google.spreadsheet.get_client(encoded_credentials)
        ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)

        assert ssheet is not None, \
            "Could not fetch %s from Google Docs." % ssheet_title

        # We allow multiple, comma-separated worksheets to be searched
        for wtitle in wsheet_title.split(','):
            wsheet = bcbio.google.spreadsheet.get_worksheet(client, ssheet, wtitle.strip())
            if not wsheet:
                logger2.warning("Could not locate %s in %s." % (wsheet_title, ssheet_title))
                continue

            # Get the rows for the project
            rows = bcbio.google.spreadsheet.get_rows_columns_with_constraint( \
                client, ssheet, wsheet, columns, {name_column: project_name})
            if len(rows) == 0:
                continue

            # Will only use the first result found to set each attribute
            for i, val in enumerate(rows[0]):
                setattr(self, attributes[i], val)

            # We have found the project data so stop iterating
            break

    @staticmethod
    def column_mapping():
        """Returns a dictionary with mapping between class attributes
        and the corresponding Google docs spreadsheet column name
        """

        return {
           "project_id": "ID",
           "project_name": "Project name",
           "queue_date": "Queue date",
           "no_samples": "No of samples",
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
