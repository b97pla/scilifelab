"""Pm deliver module"""
import os
import re
import shutil
import itertools

from cement.core import controller
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.report import sequencing_success
from scilifelab.report.rl import *
from scilifelab.report.qc import application_qc, fastq_screen, qc_cutoff
from scilifelab.report.delivery_notes import sample_status_note, project_status_note
from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection

## Main delivery controller
class DeliveryController(AbstractBaseController):
    """
    Functionality for deliveries
    """
    class Meta:
        label = 'deliver'
        description = 'Deliver data'
        arguments = [
            (['project_id'], dict(help="Project id. Standard format is 'J.Doe_00_00'", default=None, nargs="?")),
            (['flowcell_id'], dict(help="Flowcell id, formatted as AA000AAXX (i.e. without date, machine name, and run number).", default=None, nargs="?")),
            (['uppmax_project'], dict(help="Uppmax project.", default=None, nargs="?")),
            (['-i', '--interactive'], dict(help="Interactively select samples to be delivered", default=False, action="store_true")),
            (['-a', '--deliver-all-fcs'], dict(help="rsync samples from all flow cells", default=False, action="store_true")),
            ]

    @controller.expose(hide=True)
    def default(self):
        if not self._check_pargs(["project_id", "flowcell_id", "uppmax_project"]):
            return
        
## Main delivery controller
class DeliveryReportController(AbstractBaseController):
    """
    Functionality for deliveries
    """
    class Meta:
        label = 'report'
        description = 'Make delivery reports and assess qc'
        arguments = [
            (['project_id'], dict(help="Project id. Standard format is 'J.Doe_00_00'", default=None, nargs="?")),
            (['flowcell_id'], dict(help="Flowcell id, formatted as AA000AAXX (i.e. without date, machine name, and run number).", default=None, nargs="?")),
            (['-u', '--uppnex_id'], dict(help="Manually insert Uppnex project ID into the report.", default=None, action="store", type=str)),
            (['-o', '--ordered_million_reads'], dict(help="Manually insert the ordered number of read pairs (in millions)", default=None, action="store", type=str)),
            (['-r', '--customer_reference'], dict(help="Manually insert customer reference (the customer's name for the project) into reports", default=None, action="store", type=str)),
            (['--check_consistency'], dict(help="Check consistency of project sample name mapping to sample run metrics names", default=False, action="store_true")),
            (['--use_ps_map'], dict(help="Use project summary mapping in cases where no sample_run_metrics is available", default=True, action="store_false")),
            (['--use_bc_map'], dict(help="Use sample run metrics barcode mapping in cases where no sample_run_metrics is available", default=False, action="store_true")),
            (['--application'], dict(help="Set application for qc evaluation. One of '{}'".format(",".join(qc_cutoff.keys())), action="store", type=str, default=None)),
            (['--exclude_sample_ids'], dict(help="Exclude project sample ids from report generation. Provide project sample ids separated by spaces, as in '--exclude_sample_ids PS1 PS2' ", action="store", default=[], nargs="+"))
            ]

    def _process_args(self):
        pass

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Print FastQ screen output for a project/flowcell")
    def fqscreen(self):
        if not self._check_pargs(["project_id"]):
            return
        out_data = fastq_screen(**vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())


    @controller.expose(help="Print summary QC data for a flowcell/project for application QC control")
    def application_qc(self):
        if not self._check_pargs(["project_id"]):
            return
        out_data = application_qc(**vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())
            

    @controller.expose(help="Make sample status note")
    def sample_status(self):
        if not self._check_pargs(["project_id", "flowcell_id"]):
            return
        out_data = sample_status_note(**vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())

    @controller.expose(help="Make project status note")
    def project_status(self):
        if not self._check_pargs(["project_id"]):
            return
        project_status_note(**vars(self.pargs))
