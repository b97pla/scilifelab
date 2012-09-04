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

from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.utils.misc import query_yes_no, filtered_walk

## Main project controller
class ProjectController(AbstractBaseController):
    """
    Functionality for project management.
    """
    class Meta:
        label = 'project'
        description = 'Manage projects'
        arguments = [
            (['projectid'], dict(help="Scilife project id (e.g. j_doe_00_00)", default="", action="store", nargs="?")),
            (['--pbzip2'], dict(help="Use pbzip2 as compressing device", default=False, action="store_true")),
            (['--pigz'], dict(help="Use pigz as compressing device", default=False, action="store_true")),
            (['-f', '--fastq'], dict(help="Workon fastq files", default=False, action="store_true")),
            (['-p', '--pileup'], dict(help="Workon pileup files", default=False, action="store_true")),
            (['-g', '--git'], dict(help="Initialize git directory in repos and project gitdir", default=False, action="store_true")),
            (['-S', '--sampleid'], dict(help="project sample id", action="store")),
            (['-F', '--flowcellid'], dict(help="project flowcell id", action="store")),
            ]
        flowcelldir = None
        compress_opt = "-v"
        compress_prog = "gzip"
        compress_suffix = ".gz"

    ## default
    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    ## ls
    @controller.expose(help="List project folder")
    def ls(self):
        assert os.path.exists(os.path.join(self.config.get("project", "root"), self.pargs.projectid)), "no project directory %s"  % self.pargs.projectid
        if self.pargs.projectid=="":
            out = self.app.cmd.command(["ls", self.config.get("project", "root")])
        else:
            self._not_implemented("list projectid contents: only use intermediate and data directories by default" )
        if out:
            print "\n".join(self._filtered_ls(out.splitlines()))

    ## init
    @controller.expose(help="Initalize project folder")
    def init(self):
        if self.pargs.projectid=="":
            return
        self.log.info("Initalizing project %s" % self.pargs.projectid)
        ## Create directory structure
        dirs = ["%s_git" % self.pargs.projectid, "data", "intermediate"]
        gitdirs = ["config", "sbatch", "doc", "lib"] 
        [self.safe_makedir(os.path.join(self.config.get("project", "root"), self.pargs.projectid, x)) for x in dirs]
        [self.safe_makedir(os.path.join(self.config.get("project", "root"), self.pargs.projectid, dirs[0], x)) for x in gitdirs]
        ## Initialize git if repos defined and flag set
        if self.config.get("project", "repos") and self.pargs.git:
            dirs = {
                'repos':os.path.join(self.config.get("project", "repos"), "current", self.pargs.projectid),
                'gitdir':os.path.join(self.config.get("project", "root"), self.pargs.projectid, "%s_git" % self.pargs.projectid)
                    }
            self.safe_makedir(dirs['repos'])
            self.sh(["cd", dirs['repos'], "&& git init --bare"])
            self.sh(["cd", dirs['gitdir'], "&& git init && git remote add origin", dirs['repos']])

    ## utility functions
    def _assert_project(self, msg="No project defined: please supply a valid project name"):
        assert os.path.exists(os.path.join(self.config.get("project", "root"), self.pargs.projectid)), "no project directory %s"  % self.pargs.projectid
        if self.pargs.projectid=="":
            self.log.warn(msg)

    def _flowcells(self):
        self._meta.flowcelldir = os.path.join(self.config.get("project", "root"), self.pargs.projectid, "nobackup", "data")
        if not os.path.exists(self._meta.flowcelldir):
              self._meta.flowcelldir = os.path.join(self.config.get("project", "root"), self.pargs.projectid,"data")
        if not os.path.exists(self._meta.flowcelldir):
            return []
        files = os.listdir(self._meta.flowcelldir)
        return files

    ## add
    @controller.expose(help="Add boilerplate code")
    def add(self):
        self._not_implemented()

    ## clean
    @controller.expose(help="Remove files")
    def clean(self):
        self._not_implemented()

    def _set_compress_plist(self):
        ## Set pattern for compress operations
        plist = []
        if self.pargs.fastq:
            plist += [".fastq", "fastq.txt", ".fq"]
        if self.pargs.pileup:
            plist += [".pileup"]
        return plist

    def _compress(self, pattern, label="compress"):
        if self.pargs.pbzip2:
            self._meta.compress_prog = "pbzip2"
        elif self.pargs.pigz:
            self._meta.compress_prog = "pigz"

        def compress_filter(f):
            if not pattern:
                return
            return re.search(pattern, f) != None
        
        if self.pargs.input_file:
            flist = [self.pargs.input_file]
        else:
            flist = filtered_walk(os.path.join(self.config.get("project", "root"), self.pargs.projectid), compress_filter)
        
        if len(flist) > 0 and not query_yes_no("Going to {} {} files ({}...). Are you sure you want to continue?".format(label, len(flist), ",".join([os.path.basename(x) for x in flist[0:10]])), force=self.pargs.force):
            sys.exit()
        for f in flist:
            self.log.info("{}ing {}".format(label, f))
            self.app.cmd.command([self._meta.compress_prog, self._meta.compress_opt, "%s" % f], label, ignore_error=True)

    ## decompress
    @controller.expose(help="Decompress files")
    def decompress(self):
        """Decompress files"""
        self._assert_project("Not running decompress function on project root directory")
        self._meta.compress_opt = "-dv"
        if self.pargs.pbzip2:
            self._meta.compress_suffix = ".bz2"
        pattern = "|".join(["{}{}$".format(x, self._meta.compress_suffix) for x in self._set_compress_plist()])
        self._compress(pattern, label="decompress")
        
    ## compress
    @controller.expose(help="Compress files")
    def compress(self):
        self._assert_project("Not running compress function on project root directory")
        self._meta.compress_opt = "-v"
        pattern = "|".join(["{}$".format(x) for x in self._set_compress_plist()])
        self._compress(pattern)

    ## du
    @controller.expose(help="Calculate disk usage in intermediate and data directories")
    def du(self):
        self._not_implemented()

    ## deliver
    ##
    ## NOTE: this is a temporary workaround for cases where data has
    ## been removed from analysis directory
    @controller.expose(help="Deliver project data")
    def deliver(self):
        self._assert_project()
        if not self.pargs.flowcellid:
            self.log.warn("No flowcellid provided. Please provide a flowcellid from which to deliver. Available options are:\n\t{}".format("\n\t".join(self._flowcells())))
            return
        

