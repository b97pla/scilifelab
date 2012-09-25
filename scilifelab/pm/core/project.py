"""
Pm Project module
=================

Provide functionality for project management.

Commands:
^^^^^^^^^
      
       ls
         list contents
       init
         initialize a project folder
       add
         add boilerplate code
       compress
         compress files
       clean      
         remove files
       du          
         calculate disk usage
       deliver     
         deliver project data to customer

Synopsis:
---------

The following command creates a directory in the project root
named j_doe_00_00. The '-g' flag adds a git directory to the
project repos, and initializes the project subdirectory j_doe_00_00_git 
for use with git.

   pm project init j_doe_00_00 -g

FIXME: Boilerplate code can be added to the project by running

   pm project add j_doe_00_00

The boilerplate code includes makefiles, sbatch templates, and documentation
templates.

Code
====
"""
import os
import sys
import re
import yaml

from cement.core import controller, hook
from scilifelab.pm.core.controller import AbstractExtendedBaseController, AbstractBaseController
from scilifelab.pm.utils.misc import query_yes_no, filtered_walk, walk

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
    ## been removed from analysis directory
    @controller.expose(help="Transfer project data to customer. Temporary fix for cases where data has been removed from analysis directory.")
    def transfer(self):
        if not self.pargs.flowcell:
            self.log.warn("No flowcellid provided. Please provide a flowcellid from which to deliver. Available options are:\n\t{}".format("\n\t".join(self._flowcells())))
            return
        
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
            self.app.log.info("removing {}".format(f))
            self.app.cmd.safe_unlink(f)
        self.app.log.info("removing {}".format(indir))
        self.app.cmd.safe_rmdir(indir)

class BcbioRunController(AbstractBaseController):
    class Meta:
        label = 'runbcbio'
        description = 'Wrapper for bcbio analyses'
        arguments = [
            ## Also relies on the project argument from above
            (['post_process'], dict(help="post process file", action="store", default=None, nargs="?", type=str)),
            (['analysis_type'], dict(help="set analysis ", action="store", default=None, type=str, nargs="?")),
            (['--genome_build'], dict(help="genome build ", action="store", default="hg19", type=str)),
            (['--only_failed'], dict(help="only run on failed samples ", action="store_true", default=False)),
            ]
        stacked_on = 'project'

    def _sample_status(self, x):
        """Find the status of a sample.

        Look for output files: currently only look for project-summary.csv"""
        if os.path.exists(os.path.join(os.path.dirname(x), "project-summary.csv")):
            return "PASS"
        else:
            return "FAIL"

    @controller.expose(help="run automated initial analysis on samples in a project")
    def run(self):
        if not self._check_pargs(["project", "post_process", "analysis_type"]):
            return
        ## Gather sample yaml files
        pattern = "-bcbb-config.yaml$"
        flist = []
        if self.pargs.sample:
            if os.path.exists(self.pargs.sample):
                with open(self.pargs.sample) as fh:
                    flist = [x.rstrip() for x in fh.readlines()]
            else:
                pattern = "{}{}".format(self.pargs.sample, pattern)
        def bcbb_yaml_filter(f):
            return re.search(pattern, f) != None
        if not flist:
            flist = filtered_walk(os.path.join(self.app.controller._meta.project_root, self.pargs.project, "data"), bcbb_yaml_filter)
        if self.pargs.only_failed:
            status = {x:self._sample_status(x) for x in flist}
            flist = [x for x in flist if self._sample_status(x)=="FAIL"]
        if len(flist) == 0 and self.pargs.sample:
            self.app.log.info("No such sample {}".format(self.pargs.sample))
        if len(flist) > 0 and not query_yes_no("Going to start {} jobs... Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        for f in flist:
            with open(f) as fh:
                config = yaml.load(fh)
            if self.pargs.analysis_type:
                config["details"][0]["multiplex"][0]["analysis"] = self.pargs.analysis_type
                config["details"][0]["analysis"] = self.pargs.analysis_type
            if config["details"][0]["genome_build"] == 'unknown':
                config["details"][0]["genome_build"] = self.pargs.genome_build
            ## Check if files exist: if they don't, then change the suffix
            config["details"][0]["multiplex"][0]["files"].sort()
            if not os.path.exists(config["details"][0]["multiplex"][0]["files"][0]):
                if os.path.splitext(config["details"][0]["multiplex"][0]["files"][0])[1] == ".gz":
                    config["details"][0]["multiplex"][0]["files"] = [x.replace(".gz", "") for x in config["details"][0]["multiplex"][0]["files"]]
                else:
                    config["details"][0]["multiplex"][0]["files"] = ["{}.gz".format(x) for x in config["details"][0]["multiplex"][0]["files"]]
            config_file = f.replace("-bcbb-config.yaml", "-pm-bcbb-analysis-config.yaml")
            self.app.cmd.write(config_file, yaml.dump(config))
            ## Run automated_initial_analysis.py
            cur_dir = os.getcwd()
            new_dir = os.path.abspath(os.path.dirname(f))
            os.chdir(new_dir)
            self.app.cmd.command(['automated_initial_analysis.py', os.path.abspath(self.pargs.post_process), new_dir, config_file])
            os.chdir(cur_dir)

