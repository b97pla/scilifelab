"""Pm deliver module"""
import os
import re
import shutil
import itertools

from cement.core import controller
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.report import sequencing_success
from scilifelab.report.rl import *
from scilifelab.report.delivery_notes import sample_status_note, project_status_note
from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection

## QC data cutoff values
qc_cutoff = {
    'rnaseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'WG-reseq':{'PCT_PF_READS_ALIGNED':70,'PERCENT_DUPLICATION':30},
    'seqcap':{'PCT_PF_READS_ALIGNED':70,'PERCENT_ON_TARGET':60, 'PCT_TARGET_BASES_10X':90, 'PERCENT_DUPLICATION':30},
    'customcap':{'PCT_PF_READS_ALIGNED':70, 'PERCENT_DUPLICATION':30},
    'finished':{},
    }
## Mapping from genomics project list application names
application_map = {'RNA-seq (Total RNA)':'rnaseq','WG re-seq':'WG-reseq','Resequencing':'reseq', 'Exome capture':'seqcap', 'Custom':'customcap', 'Finished library':'finished' , 'Custom capture':'customcap'}
application_inv_map = {v:k for k, v in application_map.items()}

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
            (['-q', '--no-qcinfo'], dict(help="Write qcinfo to console", default=True, action="store_false")),
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
        s_con = SampleRunMetricsConnection(username=self.pargs.user, password=self.pargs.password, url=self.pargs.url)
        samples = s_con.get_samples(fc_id=self.pargs.flowcell_id, sample_prj=self.pargs.project_id)
        for s in samples:
            fqscreen_data = s.get("fastq_scr", {})
            self.app._output_data["stdout"].write(s["barcode_name"] + "\n")
            if fqscreen_data:
                header = [[x for x in v.keys()] for k, v in fqscreen_data.iteritems()]
                self.app._output_data["stdout"].write("\t\t" + "".join("{:>27}".format(x) for x in header[0]) + "\n")
                vals = ["{:>12}\t{}\n".format(k, "".join(["{:>27}".format(x) for x in v.values()])) for k, v in fqscreen_data.iteritems()]
                for v in vals:
                    self.app._output_data["stdout"].write(v)

            

    @controller.expose(help="Print summary QC data for a flowcell/project for application QC control")
    def qc(self):
        if not self._check_pargs(["project_id"]):
            return
        header = ["sample","lane","flowcell", "date",  "TOTAL_READS",
                  "MEAN_INSERT_SIZE", "GENOME_SIZE", "PERCENT_ON_TARGET",
                  "PERCENT_DUPLICATION", "PCT_TARGET_BASES_10X", "PCT_PF_READS_ALIGNED",
                  "dup_status", "status"]

        s_con = SampleRunMetricsConnection(username=self.pargs.user, password=self.pargs.password, url=self.pargs.url)
        qc_data = s_con.get_qc_data(self.pargs.project_id, self.pargs.flowcell_id)

        ## FIXME: this code is redundant. ProjectSummaryConnection is
        ## called in s_con.get_qc_data, but needed for project level
        ## info. Call get_qc_data from ProjectSummaryConnection instead?
        p_con = ProjectSummaryConnection(username=self.pargs.user, password=self.pargs.password, url=self.pargs.url)
        project = p_con.get_entry(self.pargs.project_id)
        
        if project.get("application") not in application_map.keys():
            if not self.pargs.application:
                self.app.log.warn("No such application {}. Please use the application option (available choices {})".format(app_label, ",".join(qc_cutoff.keys())))
                return
            application = self.pargs.application
        else:
            application = application_map[project.get("application")]
            
        def assess_qc(x, application):
            status = "PASS"
            dup_status = "OK"
            for k in qc_cutoff[application].keys():
                self.log.debug("assessing qc metric {}".format(k))
                if k == "PERCENT_DUPLICATION":
                    if float(x[k]) > qc_cutoff[application][k]: 
                        dup_status = "HIGH"
                else:
                    if float(x[k]) < qc_cutoff[application][k]:
                        status = "FAIL"
            genome_size = "{:.1f}G".format(int(x["GENOME_SIZE"])/1e9) if x["GENOME_SIZE"]>1e9 else "{:.1f}M".format(int(x["GENOME_SIZE"])/1e6)
            return [x["sample"],x["lane"],x["flowcell"],x["date"], 
                    "{:.2f}M".format(int(x["TOTAL_READS"])/1e6/2), "{:.1f}".format(float(x["MEAN_INSERT_SIZE"])),
                    genome_size,
                    "{:.1f}".format(float(x["PERCENT_ON_TARGET"])),
                    "{:.1f}".format(float(x["PERCENT_DUPLICATION"])),
                    "{:.1f}".format(float(x["PCT_TARGET_BASES_10X"])), 
                    "{:.1f}".format(float(x["PCT_PF_READS_ALIGNED"])), 
                    dup_status, status]
        ## Add info about qc
        self.app._output_data["stdout"].write("\n\nApplication QC\n")
        self.app._output_data["stdout"].write("==============\n")
        self.app._output_data["stdout"].write("Project\t{}\n".format(self.pargs.project_id))
        self.app._output_data["stdout"].write("Application\t{}\n".format(application_inv_map[application]))
        self.app._output_data["stdout"].write("==============\n\n")        
        self.app._output_data["stdout"].write("Application QC criteria\n")
        self.app._output_data["stdout"].write("==============\n")        
        for k,v in sorted(qc_cutoff[application].iteritems()):
            self.app._output_data["stdout"].write("{}={}\n".format(k, v))
        self.app._output_data["stdout"].write("==============\n\n")        
        self.app._output_data["stdout"].write(" ".join(header) + "\n")
        for k,v in sorted(qc_data.iteritems()):
            y = [str(x) for x in assess_qc(v, application)]
            self.app._output_data["stdout"].write("\t".join(y) + "\n")

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
