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

from cement.core import controller, hook
from pmtools.core.controller import AbstractExtendedBaseController
from pmtools.utils.misc import query_yes_no, filtered_walk

## Main project controller
class ProjectController(AbstractExtendedBaseController):
    """
    Functionality for project management.
    """
    class Meta:
        label = 'project'
        description = 'Manage projects'
        arguments = [
            (['project'], dict(help="Scilife project id (e.g. j_doe_00_00)", default="", action="store", nargs="?")),
            (['-g', '--git'], dict(help="Initialize git directory in repos and project gitdir", default=False, action="store_true")),
            (['-S', '--sample'], dict(help="project sample id", action="store")),
            (['-F', '--flowcell'], dict(help="project flowcell id", action="store")),
            (['--finished'], dict(help="include finished project listing", action="store_true", default=False)),
            ]
        flowcelldir = None

    def _process_args(self):
        # setup project search space
        if self.app.pargs.finished:
            self._meta.project_root = self.app.config.get("project", "finished")
        else:
            self._meta.project_root = self.app.config.get("project", "root")
        # Set root path for parent class
        self._meta.root_path = self._meta.project_root
        assert os.path.exists(self._meta.project_root), "No such directory {}; check your project config".format(self._meta.project_root)
        ## Set path_id for parent class
        if self.pargs.project:
            self._meta.path_id = self.pargs.project
        super(ProjectController, self)._process_args()

    ## utility functions
    def _assert_project_id(self, msg):
        if not self.pargs.project:
            self.log.warn(msg)
            return False
        return True

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
    @controller.expose(help="Transfer project data")
    def transfer(self):
        if not self.pargs.flowcell:
            self.log.warn("No flowcellid provided. Please provide a flowcellid from which to deliver. Available options are:\n\t{}".format("\n\t".join(self._flowcells())))
            return
        

