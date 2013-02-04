"""Pm deliver module"""
import os
import re
import shutil
import itertools
from cement.core import controller
from scilifelab.pm.core.controller import AbstractBaseController, AbstractExtendedBaseController
from scilifelab.report import sequencing_success
from scilifelab.report.rl import *
from scilifelab.report.qc import application_qc, fastq_screen, QC_CUTOFF
from scilifelab.bcbio.run import find_samples
from scilifelab.report.delivery_notes import sample_status_note, project_status_note
from scilifelab.report.best_practice import best_practice_note, SEQCAP_KITS
from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, get_scilife_to_customer_name
from scilifelab.utils.misc import query_yes_no, filtered_walk

BCBIO_EXCLUDE_DIRS = ['realign-split', 'variants-split', 'tmp', 'tx', 'fastqc', 'fastq_screen', 'alignments', 'nophix']

## Main delivery controller
class DeliveryController(AbstractBaseController):
    """
    Functionality for deliveries
    """
    class Meta:
        label = 'deliver'
        description = 'Deliver data'
        root_path = None
        path_id = None
        project_root = None
        arguments = [
            (['project'], dict(help="Project name, formatted as 'J.Doe_00_00'", default=None, nargs="?")),
            (['uppmax_project'], dict(help="Uppmax project.", default=None, nargs="?")),
            (['-i', '--interactive'], dict(help="Interactively select samples to be delivered", default=False, action="store_true")),
            (['-a', '--deliver_all_fcs'], dict(help="rsync samples from all flow cells", default=False, action="store_true")),
            (['-S', '--sample'], dict(help="Project sample id. If sample is a file, read file and use sample names within it. Sample names can also be given as full paths to bcbb-config.yaml configuration file.", action="store", default=None, type=str)),
            (['--no_bam'], dict(help="Don't include bam files in delivery", action="store", default=False)),
            (['--no_vcf'], dict(help="Don't include vcf files in delivery", action="store", default=False)),
            ]

    def _setup(self, base_app):
        super(DeliveryController, self)._setup(base_app)
        group = base_app.args.add_argument_group('delivery', 'Options affecting data delivery.')
        group.add_argument('--move', help="Transfer file with move", default=False, action="store_true")
        group.add_argument('--copy', help="Transfer file with copy", default=False, action="store_true")
        group.add_argument('--rsync', help="Transfer file with rsync (default)", default=True, action="store_true")
        group.add_argument('--intermediate', help="Work on intermediate data", default=False, action="store_true")
        group.add_argument('--data', help="Work on data folder", default=False, action="store_true")

    def _process_args(self):
        # NB: duplicate of project.ProjectController._process_args
        # setup project search space
        self._meta.project_root = self.app.config.get("project", "root")
        # Set root path for parent class
        self._meta.root_path = self._meta.project_root
        assert os.path.exists(self._meta.project_root), "No such directory {}; check your project config".format(self._meta.project_root)
        if self.pargs.project:
            self._meta.path_id = self.pargs.project
            # Add intermediate or data
            if self.app.pargs.intermediate:
                if os.path.exists(os.path.join(self._meta.project_root, self._meta.path_id, "nobackup")):
                    self._meta.path_id = os.path.join(self._meta.path_id, "nobackup", "intermediate")
                else:
                    self._meta.path_id = os.path.join(self._meta.path_id, "intermediate")
            if self.app.pargs.data and not self.app.pargs.intermediate:
                if os.path.exists(os.path.join(self._meta.project_root, self._meta.path_id, "nobackup")):
                    self._meta.path_id = os.path.join(self._meta.path_id, "nobackup", "data")
                else:
                    self._meta.path_id = os.path.join(self._meta.path_id, "data")
        # Setup transfer options
        if self.pargs.move:
            self.pargs.rsync = False
        elif self.pargs.copy:
            self.pargs.rsync = False
        super(DeliveryController, self)._process_args()

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Deliver best practice results")
    def best_practice(self):
        if not self._check_pargs(["project", "uppmax_project"]):
            return
        outpath = os.path.normpath(os.path.abspath(self.pargs.uppmax_project))
        if not os.path.exists(outpath):
            self.app.cmd.safe_makedir(outpath)
        kw = vars(self.pargs)
        basedir = os.path.abspath(os.path.join(self._meta.root_path, self._meta.path_id))
        flist = find_samples(basedir, **vars(self.pargs))
        if not len(flist) > 0:
            self.log.info("No samples/sample configuration files found")
            return
        def filter_fn(f):
            if not pattern:
                return
            return re.search(pattern, f) != None
        # Setup pattern
        plist = [".*.yaml$", ".*.metrics$"]
        if not self.pargs.no_bam:
            plist.append(".*.bam$")
        if not self.pargs.no_vcf:
            plist.append(".*.vcf$")
        pattern = "|".join(plist)
        for f in flist:
            path = os.path.dirname(f)
            sources = filtered_walk(path, filter_fn=filter_fn, exclude_dirs=BCBIO_EXCLUDE_DIRS)
            targets = [src.replace(basedir, outpath) for src in sources]
            self._transfer_files(sources, targets)

    def _transfer_files(self, sources, targets):
        for src, tgt in zip(sources, targets):
            if not os.path.exists(os.path.dirname(tgt)):
                self.app.cmd.safe_makedir(os.path.dirname(tgt))
            self.app.cmd.transfer_file(src, tgt)
            
## Main delivery controller
class DeliveryReportController(AbstractBaseController):
    """
    Functionality for deliveries
    """
    class Meta:
        label = 'report'
        description = 'Make delivery reports and assess qc'

    def _setup(self, app):
        group = app.args.add_argument_group('Reporting options', 'Options that affect report output')
        group.add_argument('project_name', help="Project name. Standard format is 'J.Doe_00_00'", default=None, nargs="?")
        group.add_argument('flowcell', help="Flowcell id, formatted as AA000AAXX (i.e. without date, machine name, and run number).", default=None, nargs="?")
        group.add_argument('-u', '--uppnex_id', help="Manually insert Uppnex project ID into the report.", default=None, action="store", type=str)
        group.add_argument('-o', '--ordered_million_reads', help="Manually insert the ordered number of read pairs (in millions), either as a string to set all samples, or a JSON string or JSON file to set at a sample level.", default=None, action="store", type=str)
        group.add_argument('-r', '--customer_reference', help="Manually insert customer reference (the customer's name for the project) into reports", default=None, action="store", type=str)
        group.add_argument('--application', help="Set application for qc evaluation. One of '{}'".format(",".join(QC_CUTOFF.keys())), action="store", type=str, default=None)
        group.add_argument('--exclude_sample_ids', help="Exclude project sample ids from report generation. Input is either a string or a JSON file with a key:value mapping, as in '--exclude_sample_ids \"{'PS1':[], 'PS2':['AACCGG']}\"'. The values consist of a list of barcodes; if the list is empty, exclude the entire sample.", action="store", default={})
        group.add_argument('--bc_count', help="Manually set barcode counts in *millions of reads*. Input is either a string or a JSON file with a key:value mapping, as in '--bc_count \"{'Sample1':100, 'Sample2':200}\"'.", action="store", default={})
        group.add_argument('--sample_aliases', help="Provide sample aliases for cases where project summary has multiple names for a sample. Input is either a string or a JSON file with a key:value mapping, for example '--sample_aliases \"{'sample1':['alias1_1', 'alias1_2'], 'sample2':['alias2_1']}\", where the value is a list of aliases. The key will be used as 'base' information, possibly updated by information from the alias entry.", action="store", default={})
        group.add_argument('--project_alias', help="Provide project aliases for cases where project summary has multiple names for a project. Input is a comma-separated list of names enclosed by brackets, for example '--project_alias \"['alias1']\"", action="store", default=None)
        group.add_argument('--phix', help="Provide phix error rate for new illumina flowcells where phix error rate is missing. Input is either a string, or a dictionary/JSON file with a lane:error mapping, for example '--phix \"{1:0.3, 2:0.4}\".", action="store", default=None)
        group.add_argument('--sphinx', help="Generate editable sphinx template. Installs conf.py and Makefile for subsequent report generation.", action="store", default=None, type=float)
        group.add_argument('--project_id', help="Project identifier, formatted as 'P###'.",  action="store", default=None, type=str)
        group.add_argument('--include_all_samples', help="Include all samples in project status report. Default is to only use the latest library prep.",  action="store_true", default=False)
        super(DeliveryReportController, self)._setup(app)

    def _process_args(self):
        pass

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Print FastQ screen output for a project/flowcell")
    def fqscreen(self):
        if not self._check_pargs(["project_name"]):
            return
        out_data = fastq_screen(**vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())


    @controller.expose(help="Print summary QC data for a flowcell/project for application QC control")
    def application_qc(self):
        if not self._check_pargs(["project_name"]):
            return
        out_data = application_qc(**vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())
            

    @controller.expose(help="Make sample status note")
    def sample_status(self):
        if not self._check_pargs(["project_name", "flowcell"]):
            return
        kw = vars(self.pargs)
        kw.update({"samplesdb":self.app.config.get("db", "samples"), "flowcelldb":self.app.config.get("db", "flowcells"), "projectdb":self.app.config.get("db", "projects")})
        out_data = sample_status_note(**kw)
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())
        self.app._output_data['debug'].write(out_data['debug'].getvalue())

    @controller.expose(help="Make project status note")
    def project_status(self):
        if not self._check_pargs(["project_name"]):
            return
        kw = vars(self.pargs)
        kw.update({"samplesdb":self.app.config.get("db", "samples"), "flowcelldb":self.app.config.get("db", "flowcells"), "projectdb":self.app.config.get("db", "projects")})
        out_data = project_status_note(**kw)
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())
        self.app._output_data['debug'].write(out_data['debug'].getvalue())
        
    @controller.expose(help="Make best practice reports")
    def best_practice(self):
        self.log.info("Until best practice results are stored in statusDB, best practice reports are generated via the 'pm project bpreport' subcommand. This requires that best practice analyses have been run in the project folder.")
        return

class BestPracticeReportController(AbstractBaseController):
    class Meta:
        label = 'bpreport'
        description = 'Functions for generating best practice reports'

    def _setup(self, app):
        group = app.args.add_argument_group('Best practice report group', 'Options for bpreport')
        group.add_argument('--application', help="Set application for best practice application note." , action="store", type=str, default="seqcap")
        group.add_argument('--capture_kit', help="Set capture for seqcap application note. Either a key, one of '{}', or free text description.".format(",".join(SEQCAP_KITS.keys())), action="store", type=str, default="agilent_v4")
        group.add_argument('--no_statusdb', help="Don't statusdb to convert scilife names to customer names.", action="store_true", default=False)
        group.add_argument('--statusdb_project_name', help="Project name in statusdb.", action="store", default=None)
        super(BestPracticeReportController, self)._setup(app)

    @controller.expose(help="Make best practice reports")
    def bpreport(self):
        if not self._check_pargs(["project"]):
            return
        kw = vars(self.pargs)
        basedir = os.path.abspath(os.path.join(self.app.controller._meta.root_path, self.app.controller._meta.path_id))
        flist = find_samples(basedir, **vars(self.pargs))
        if not len(flist) > 0:
            self.log.info("No samples/sample configuration files found")
            return
        if self.pargs.no_statusdb:
            sample_name_map = None
        else:
            if not self._check_pargs(["statusdb_project_name"]):
                return
            p_con = ProjectSummaryConnection(dbname=self.app.config.get("db", "projects"), **vars(self.app.pargs))
            s_con = SampleRunMetricsConnection(dbname=self.app.config.get("db", "samples"), **vars(self.app.pargs))
            sample_name_map = get_scilife_to_customer_name(self.pargs.statusdb_project_name, p_con, s_con)
        kw.update(project_name=self.pargs.project, flist=flist, basedir=basedir, sample_name_map=sample_name_map)
        out_data = best_practice_note(**kw)
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())
        self.app._output_data['debug'].write(out_data['debug'].getvalue())


