"""Create reports on google docs
"""
# import copy
import logbook
import datetime
import yaml

from bcbio.google import _from_unicode
from bcbio.google import get_credentials
from bcbio.google import bc_metrics
from bcbio.google import qc_metrics
from bcbio.google import document as g_document
from bcbio.google import spreadsheet as g_spreadsheet
from bcbio.pipeline.qcsummary import RTAQCMetrics
from bcbio.pipeline.flowcell import Flowcell
from bcbio.log import create_log_handler
from bcbio.log import logger2 as log

from bcbio.distributed import messaging
import os


def queue_report(fc_date, fc_name, run_info_yaml, dirs, config, config_file):
    if "gdocs_upload" not in config:
        return False

    runner = messaging.runner("bcbio.distributed.google_tasks", \
        {"work": os.getcwd(), \
         "config": os.path.dirname(config_file)}, \
        config, \
        config_file, \
        wait=False)

    runner("create_report_on_gdocs", [[fc_date, fc_name, run_info_yaml, dirs, config]])

    return True


def create_report_on_gdocs(fc_date, fc_name, run_info_yaml, dirs, config):
    """Create reports on gdocs containing both demultiplexed read counts and QC data.
    """
    success = True
    try:
        # Inject the fc_date and fc_name in the email subject
        def record_processor(record):
            return record.extra.__setitem__('run', "%s_%s" % (fc_date, fc_name))

        # Parse the run_info.yaml file
        log.debug("Loading this run_info: {}".format(run_info_yaml))
        with open(run_info_yaml, "r") as fh:
            run_info = yaml.load(fh)

        # Get the gdocs account credentials
        encoded_credentials = get_credentials(config)
        if not encoded_credentials:
            log.warn("Could not find Google Docs account credentials in configuration. \
                      No sequencing report was written")
            return False

        # Get the required parameters from the post_process.yaml configuration file
        gdocs = config.get("gdocs_upload", None)

        # Add email notification
        email = gdocs.get("gdocs_email_notification", None)
        smtp_host = config.get("smtp_host", "")
        smtp_port = config.get("smtp_port", "")
        log_handler = create_log_handler({'email': email, \
                                          'smtp_host': smtp_host, \
                                          'smtp_port': smtp_port}, True)

    except Exception as e:
        success = False
        log.warn("Encountered exception when writing sequencing report to Google Docs: %s" % e)

    with log_handler.applicationbound(), logbook.Processor(record_processor):
        try:
            log.info("Started creating sequencing report on Google docs for %s_%s on %s" \
                % (fc_date, fc_name, datetime.datetime.now().isoformat()))

            # Get a flowcell object
            fc = Flowcell(fc_name, fc_date, run_info, dirs.get("work", None))

            # Get the GDocs demultiplex result file title
            gdocs_dmplx_spreadsheet = gdocs.get("gdocs_dmplx_file", None)
            # Get the GDocs QC file title
            gdocs_qc_spreadsheet = gdocs.get("gdocs_qc_file", None)

            # FIXME: Make the bc stuff use the Flowcell module
            if gdocs_dmplx_spreadsheet is not None:
                # Upload the data
                bc_metrics.write_run_report_to_gdocs(fc, fc_date, \
                    fc_name, gdocs_dmplx_spreadsheet, encoded_credentials, append=True)
            else:
                log.warn("Could not find Google Docs demultiplex results file \
                    title in configuration. No demultiplex counts were \
                    written to Google Docs for %s_%s" % (fc_date, fc_name))

            # Parse the QC metrics
            try:
                qc = RTAQCMetrics(dirs.get("flowcell", None))
            except:
                qc = None

            if gdocs_qc_spreadsheet is not None and qc is not None:
                qc_metrics.write_run_report_to_gdocs(fc, qc, gdocs_qc_spreadsheet, encoded_credentials)
            else:
                log.warn("Could not find Google Docs QC file title in configuration. " \
                         "No QC data were written to Google Docs " \
                         "for %s_%s".format(fc_date, fc_name))

            # Get the projects parent folder
            projects_folder = gdocs.get("gdocs_projects_folder", None)

            # Write the bc project summary report
            if projects_folder is not None:
                create_project_report_on_gdocs(fc, qc, \
                    encoded_credentials, projects_folder)

        except Exception as e:
            success = False
            log.warn("Encountered exception when writing sequencing report " \
                     "to Google Docs: {}".format(e))

        if success:
            log.info("Sequencing report successfully created on Google " \
                     "docs for {}_{} on {}".format(fc_date, fc_name, datetime.datetime.now().isoformat()))
        else:
            log.warn("Encountered exception when writing sequencing " \
                     "report for %s_%s to Google docs on %s" \
                     % (fc_date, fc_name, datetime.datetime.now().isoformat()))

    return success


def create_project_report_on_gdocs(fc, qc, encoded_credentials, gdocs_folder):
    """Upload the sample read distribution for a project to google docs.
    """
    success = True
    # Create a client class which will make HTTP requests with Google Docs server.
    client = g_spreadsheet.get_client(encoded_credentials)
    doc_client = g_document.get_client(encoded_credentials)

    # Get a reference to the parent folder
    parent_folder = g_document.get_folder(doc_client, gdocs_folder)

    if not parent_folder:
        parent_folder = g_document.add_folder(doc_client, gdocs_folder)
        log.info("Folder {!r} created".format(gdocs_folder))

    parent_folder_title = _from_unicode(parent_folder.title.text)

    # Loop over the projects
    for project_name in fc.get_project_names():

        # Get a flowcell object containing just the data for the project
        project_fc = fc.prune_to_project(project_name)

        ssheet_title = project_name + "_sequencing_results"
        ssheet = g_spreadsheet.get_spreadsheet(client, ssheet_title)
        if not ssheet:
            ssheet = g_document.add_spreadsheet(doc_client, ssheet_title)
            g_document.move_to_folder(doc_client, ssheet, parent_folder)
            ssheet = g_spreadsheet.get_spreadsheet(client, ssheet_title)
            log.info("Spreadsheet {!r} created in " \
                     "folder {!r}".format(_from_unicode(ssheet.title.text), \
                                          parent_folder_title))

        success &= bc_metrics._write_project_report_to_gdocs(client, ssheet, project_fc)
        success &= bc_metrics.write_project_report_summary_to_gdocs(client, ssheet)
        success &= qc_metrics.write_run_report_to_gdocs(project_fc, \
                                                        qc, \
                                                        ssheet_title, \
                                                        encoded_credentials)

        log.info("Sequencing results report written to spreadsheet '%s'" \
            % _from_unicode(ssheet.title.text))

    return success
