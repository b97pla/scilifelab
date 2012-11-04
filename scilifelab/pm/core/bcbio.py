"""pm bcbio extension"""
import os
import re

from cement.core import backend, controller, handler, hook
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.misc import query_yes_no, filtered_walk
from scilifelab.bcbio.run import find_samples, setup_sample, remove_files, run_bcbb_command, setup_merged_samples

class BcbioRunController(AbstractBaseController):
    class Meta:
        label = 'bcbio'
        description = 'Wrapper for bcbio analyses'

    def _setup(self, app):
        group = app.args.add_argument_group('Bcbio argument group', 'Options for bcbio')
        group.add_argument('-M', '--merge_analysis', action='store_true',
                       help='Run merge analysis. If a sample has data from more than one run, a \'total\' directory will be setup.', default=False)
        group.add_argument('--post_process', help="post process file. Setting this will override sample-specific post process files. Currently not implemented.", action="store", default=None, nargs="?", type=str)
        group.add_argument('--genome_build', help="genome build ", action="store", default=None, type=str)
        group.add_argument('--only_failed', help="only run on failed samples ", action="store_true", default=False)
        group.add_argument('--amplicon', help="amplicon-based analyses (e.g. HaloPlex), which means mark_duplicates is set to false", action="store_true", default=False)
        group.add_argument('--targets', help="sequence capture target file", action="store", default=None)
        group.add_argument('--baits', help="sequence capture baits file", action="store", default=None)
        group.add_argument('--distributed', help="run distributed, changing 'num_cores' in  post_process to 'messaging': calls automated_initial_analysis.py", action="store_true", default=False)
        group.add_argument('--num_cores', help="num_cores value; default 8", action="store", default=8, type=int)
        group.add_argument('--only_setup', help="only perform setup", action="store_true", default=False)
        group.add_argument('--restart', help="restart analysis", action="store_true", default=False)
        group.add_argument('--analysis', help="set analysis type in bcbb config file", action="store", default=None, type=str)
        group.add_argument('--snpEff', help="set snpEff program in post process", action="store", default=None, type=str)
        group.add_argument('--no_merged', help="dot't do merged sample analysis for samples with multiple sample runs", action="store_true", default=False)
        super(BcbioRunController, self)._setup(app)

    def _sample_status(self, x):
        """Find the status of a sample.

        Look for output files: currently only look for project-summary.csv"""
        if os.path.exists(os.path.join(os.path.dirname(x), "project-summary.csv")):
            return "PASS"
        else:
            return "FAIL"

    @controller.expose(help="Run bcbb pipeline")
    def run(self):
        if not self._check_pargs(["project"]):
            return
        flist = find_samples(os.path.abspath(os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id)), **vars(self.pargs))
        if not self.pargs.no_merged:
            ##  Setup merged samples and append to flist if new list longer
            flist = setup_merged_samples(flist, **vars(self.pargs))
        if len(flist) > 0 and not query_yes_no("Going to start {} jobs... Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        orig_dir = os.path.abspath(os.getcwd())
        for run_info in flist:
            os.chdir(os.path.abspath(os.path.dirname(run_info)))
            setup_sample(run_info, **vars(self.pargs))
            os.chdir(orig_dir)
        if self.pargs.only_setup:
            return
        if self.pargs.only_failed:
            status = {x:self._sample_status(x) for x in flist}
            flist = [x for x in flist if self._sample_status(x)=="FAIL"]
        ## Here process files again, removing if requested, and running the pipeline
        for run_info in flist:
            self.app.log.info("Running analysis defined by config file {}".format(run_info))
            os.chdir(os.path.abspath(os.path.dirname(run_info)))
            if self.app.cmd.monitor(work_dir=os.path.dirname(run_info)):
                self.app.log.warn("Not running job")
                continue
            if self.pargs.restart:
                self.app.log.info("Removing old analysis files in {}".format(os.path.dirname(run_info)))
                remove_files(run_info, **vars(self.pargs))
            (cl, platform_args) = run_bcbb_command(run_info, **vars(self.pargs))
            self.app.cmd.command(cl, **{'platform_args':platform_args, 'saveJobId':True})
            os.chdir(orig_dir)
