"""
Pm Project module

Provide functionality for project management.

Commands:
       ls          list contents
       init        initialize a project folder
       add         add boilerplate code
       compress    compress files
       clean       remove files
       du          calculate disk usage
       deliver     deliver project data to customer

Synopsis:

The following command creates a directory in the project root
named j_doe_00_00. The '-g' flag adds a git directory to the
project repos, and initializes the project subdirectory j_doe_00_00_git 
for use with git.

   pm project init j_doe_00_00 -g

FIXME: Boilerplate code can be added to the project by running

   pm project add j_doe_00_00

The boilerplate code includes makefiles, sbatch templates, and documentation
templates.
"""
import os
import sys
import re

from cement.core import controller
from pmtools import AbstractBaseController

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
            (['-s', '--sbatch'], dict(help="Submit jobs to slurm", default=False, action="store_true")),
            (['-f', '--fastq'], dict(help="Workon fastq files", default=False, action="store_true")),
            (['-p', '--pileup'], dict(help="Workon pileup files", default=False, action="store_true")),
            (['-A', '--uppmax-project'], dict(help="uppmax project id for use with sbatch", action="store")),
            (['-t', '--sbatch-time'], dict(help="sbatch time limit", default="00:10:00", action="store")),
            (['-N', '--node'], dict(help="run node job", default=False, action="store_true")),
            (['-g', '--git'], dict(help="Initialize git directory in repos and project gitdir", default=False, action="store_true")),
            (['-S', '--sampleid'], dict(help="project sample id", action="store")),
            (['-F', '--flowcellid'], dict(help="project flowcell id", action="store")),
            ]

    ## default
    @controller.expose(hide=True)
    def default(self):
        print __doc__

    ## ls
    @controller.expose(help="List project folder")
    def ls(self):
        assert os.path.exists(os.path.join(self.config.get("project", "root"), self.pargs.projectid)), "no project directory %s"  % self.pargs.projectid
        if self.pargs.projectid=="":
            out = self.sh(["ls", self.config.get("project", "root")])
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
            print self.log._meta.namespace
            print "in assert"

            print dir(self.log.app)
            self.log.info(msg)
            self.log.warn(msg)
            sys.exit()
            sys.exit()

    def _flowcells(self):
        files = os.listdir(os.path.join(self.config.get("project", "root")), "data")
        return files

    ## add
    @controller.expose(help="Add boilerplate code")
    def add(self):
        self._not_implemented()

    ## clean
    @controller.expose(help="Remove files")
    def clean(self):
        self._not_implemented()

    ## compress
    @controller.expose(help="Compress files")
    def compress(self):
        self._assert_project("Not running compress function on project root directory")
        ## Set pattern for compress operations
        plist = []
        if self.pargs.fastq:
            plist += [".fastq$", ".fastq.txt$", ".fq$"]
        if self.pargs.pileup:
            plist += [".pileup$"]
        pattern = "|".join(plist)

        def compress_filter(f):
            return re.search(pattern, f) != None
        flist = []
        for root, dirs, files in os.walk(os.path.join(self.config.get("project", "root"), self.pargs.projectid)):
            flist = flist + [os.path.join(root, x) for x in filter(compress_filter, files)]
        ##self.log.info("Going to compress %s files. Are you sure you want to continue?" % len(flist))
        if not self.query_yes_no("Going to compress %s files. Are you sure you want to continue?" % len(flist)):
            sys.exit()
        for f in flist:
            print "compressing %s" %f
            self.drmaa(["gzip",  "-v",  "%s" % f], "compress")

    ## du
    @controller.expose(help="Calculate disk usage in intermediate and data directories")
    def du(self):
        self._not_implemented()


    ## deliver
    ##
    ## NOTE: this is a temporary workaround for cases where data has
    ## been removed from analysis directory
    @controller.expose(help="Deliver project data")
    def delivery(self):
        self._assert_project()
        ## FIXME: list a flowcell and deliver from there
        print self._flowcells()

