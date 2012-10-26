"""Pm Project module"""
import os
import sys
import re
import yaml

from cement.core import controller, hook
from scilifelab.pm.core.controller import AbstractExtendedBaseController, AbstractBaseController
from scilifelab.utils.misc import query_yes_no, filtered_walk, walk
from scilifelab.pm.lib.clean import purge_alignments
from scilifelab.bcbio.run import find_samples, setup_sample, remove_files, run_bcbb_command


## Main project controller
class ProjectController(AbstractExtendedBaseController):
    """
    Functionality for project management.
    """
    class Meta:
        label = 'project'
        description = 'Manage projects'
        arguments = [
            (['project'], dict(help="Scilife project id (e.g. j_doe_00_00)", default=None, action="store", nargs="?")),
            (['-g', '--git'], dict(help="Initialize git directory in repos and project gitdir", default=False, action="store_true")),
            (['-S', '--sample'], dict(help="project sample id", action="store", default=None, type=str)),
            (['-F', '--flowcell'], dict(help="project flowcell id", action="store", default=None, type=str)),
            (['--finished'], dict(help="include finished project listing", action="store_true", default=False)),
            (['--intermediate'], dict(help="Work on intermediate data", default=False, action="store_true")),
            (['--data'], dict(help="Work on data folder", default=False, action="store_true")),
            (['--minfilesize'], dict(help="Min file size to keep (in bytes). Default 2000.", default=2000, action="store", type=int)),
            ]
        flowcelldir = None

    ## Remember: need to do argument processing here also for stacked controllers
    ## FIX ME: _process_args should be called in the stacked controller
    def _process_args(self):
        # setup project search space
        if self.app.pargs.finished:
            self._meta.project_root = self.app.config.get("project", "finished")
        else:
            self._meta.project_root = self.app.config.get("project", "root")

        # If rm function set intermediate
        if self.command == "rm":
            self.pargs.intermediate = True

        # Set root path for parent class
        self._meta.root_path = self._meta.project_root
        assert os.path.exists(self._meta.project_root), "No such directory {}; check your project config".format(self._meta.project_root)
        ## Set path_id for parent class
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
        super(ProjectController, self)._process_args()
 
    ## init
    @controller.expose(help="Initalize project folder")
    def init(self):
        if self.pargs.project=="":
            return
        self.log.info("Initalizing project %s" % self.pargs.project)
        ## Create directory structure
        dirs = ["%s_git" % self.pargs.project, "data", "intermediate"]
        gitdirs = ["config", "sbatch", "doc", "lib"] 
        [self.safe_makedir(os.path.join(self._meta.project_root, self.pargs.project, x)) for x in dirs]
        [self.safe_makedir(os.path.join(self._meta.project_root, self.pargs.project, dirs[0], x)) for x in gitdirs]
        ## Initialize git if repos defined and flag set
        if self.config.get("project", "repos") and self.pargs.git:
            dirs = {
                'repos':os.path.join(self.config.get("project", "repos"), "current", self.pargs.project),
                'gitdir':os.path.join(self._meta.project_root, self.pargs.project, "%s_git" % self.pargs.project)
                    }
            self.safe_makedir(dirs['repos'])
            self.sh(["cd", dirs['repos'], "&& git init --bare"])
            self.sh(["cd", dirs['gitdir'], "&& git init && git remote add origin", dirs['repos']])

    def _flowcells(self):
        files = []
        if self.pargs.project:
            self._meta.flowcelldir = os.path.join(self._meta.project_root, self.pargs.project, "nobackup", "data")
            if not os.path.exists(self._meta.flowcelldir):
                self._meta.flowcelldir = os.path.join(self._meta.project_root, self.pargs.project,"data")
            if not os.path.exists(self._meta.flowcelldir):
                return []
            files = os.listdir(self._meta.flowcelldir)
        return files

    ## add
    @controller.expose(help="Add boilerplate code")
    def add(self):
        self._not_implemented()
        
    ## NOTE: this is a temporary workaround for cases where data has
    ## been removed from production directory
    @controller.expose(help="Transfer project data to customer. Temporary fix for cases where data has been removed from production directory.")
    def transfer(self):
        if not self.pargs.flowcell:
            self.log.warn("No flowcellid provided. Please provide a flowcellid from which to deliver. Available options are:\n\t{}".format("\n\t".join(self._flowcells())))
            return
    
    ## purge_alignments
    @controller.expose(help="purge alignments in project folders")
    def purge(self):
        """Cleanup sam and bam files. In some cases, sam files
        persist. If the corresponding bam file exists, replace the sam
        file contents with a message that the file has been removed to
        save space.
        """
        if not self._check_pargs(["project"]):
            return
        if self.app.pargs.sam:
            purge_alignments(path=os.path.join(self._meta.root_path, self._meta.path_id), dry_run=self.app.pargs.dry_run, force=self.app.pargs.force, fsize=self.app.pargs.minfilesize)
        else:
            purge_alignments(path=os.path.join(self._meta.root_path, self._meta.path_id), dry_run=self.app.pargs.dry_run, force=self.app.pargs.force, ftype="bam", fsize=self.app.pargs.minfilesize)

class ProjectRmController(AbstractBaseController):
    class Meta:
        label = 'projectrm'
        description = 'Functionality for removing analyses from intermediate folders'
        arguments = [
            (['analysis_id'], dict(help="analysis name in intermediate", action="store", default=None, nargs="?", type=str)),
            ]
        stacked_on = 'project'

    @controller.expose(help="Remove analyses from project intermediate subfolder")
    def rm(self):
        if not self._check_pargs(["project",  "analysis_id"]):
            return
        indir = os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id, self.pargs.analysis_id)
        assert os.path.exists(indir), "No such analysis {} for project {}".format(self.pargs.analysis_id, self.pargs.project)
        try:
            flist = walk(indir)
        except IOError as e:
            self.app.log.warn(str(e))
            raise e
        if len(flist) > 0 and not query_yes_no("Going to remove all contents ({} files) of analysis {} for project {}... Are you sure you want to continue?".format(len(flist), self.pargs.analysis_id, self.pargs.project), force=self.pargs.force):
            return
        for f in flist:
            self.app.cmd.safe_unlink(f)
        self.app.log.info("removing {}".format(indir))
        self.app.cmd.safe_rmdir(indir)


class BcbioRunController(AbstractBaseController):
    class Meta:
        label = 'runbcbio'
        description = 'Wrapper for bcbio analyses'
        arguments = [
            ## Also relies on the project argument from above
            (['--post_process'], dict(help="post process file. Setting this will override sample-specific post process files. Currently not implemented.", action="store", default=None, nargs="?", type=str)),
            (['--genome_build'], dict(help="genome build ", action="store", default="hg19", type=str)),
            (['--only_failed'], dict(help="only run on failed samples ", action="store_true", default=False)),
            (['--amplicon'], dict(help="amplicon-based analyses (e.g. HaloPlex), which means mark_duplicates is set to false", action="store_true", default=False)),
            (['--targets'], dict(help="sequence capture target file", action="store", default=None)),
            (['--baits'], dict(help="sequence capture baits file", action="store", default=None)),
            (['--distributed'], dict(help="run distributed, changing 'num_cores' in  post_process to 'messaging': calls automated_initial_analysis.py", action="store_true", default=False)),
            (['--num_cores'], dict(help="num_cores value; default 8", action="store", default=8, type=int)),
            (['--only_setup'], dict(help="only perform setup", action="store_true", default=False)),
            (['--restart'], dict(help="restart analysis", action="store_true", default=False)),
            (['--no_only_run'], dict(help="run_bcbb parameter: don't setup", action="store_true", default=False)),
            (['--google_report'], dict(help="make a google report (default False)", action="store_false", default=True)),
            (['--analysis_type'], dict(help="set analysis type in bcbb config file", action="store", default="Align_standard_seqcap", type=str)),
            (['--email'], dict(help="set user email address", action="store", default=None, type=str)),
            ]
        stacked_on = 'project'

    def _sample_status(self, x):
        """Find the status of a sample.

        Look for output files: currently only look for project-summary.csv"""
        if os.path.exists(os.path.join(os.path.dirname(x), "project-summary.csv")):
            return "PASS"
        else:
            return "FAIL"

    ## FIXME: for now this is an exact copy of production.run. Since
    ## this is a function related to bcbio, it should be put in an
    ## extension ext_bcbio.py. Problem is: how can a controller be
    ## stacked on two other controllers, namely production and
    ## project?
    @controller.expose(help="Run bcbb pipeline")
    def run(self):
        if not self._check_pargs(["project"]):
            return
        flist = find_samples(os.path.abspath(os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id)), **vars(self.pargs))
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
