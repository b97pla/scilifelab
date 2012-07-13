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
            ]

    @controller.expose(hide=True)
    def default(self):
        __doc__

    @controller.expose(help="List project folder")
    def ls(self):
        assert os.path.exists(os.path.join(self.config.get("project", "root"), self.pargs.projectid)), "no project directory %s"  % self.pargs.projectid
        if self.pargs.projectid=="":
            out = self.sh(["ls", self.config.get("project", "root")])
        else:
            self._not_implemented("list projectid contents: only use intermediate and data directories by default" )
        if out:
            print "\n".join(self._filtered_ls(out.splitlines()))

    @controller.expose(help="Initalize project folder")
    def init(self):
        self._not_implemented()

    @controller.expose(help="Add boilerplate code")
    def add(self):
        self._not_implemented()

    @controller.expose(help="Remove files")
    def clean(self):
        self._not_implemented()

    @controller.expose(help="Compress files")
    def compress(self):
        assert os.path.exists(os.path.join(self.config.get("project", "root"), self.pargs.projectid)), "no project directory %s"  % self.pargs.projectid
        if self.pargs.projectid=="":
            self.log.warn("Not running compress function on project root directory")
            sys.exit()

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
            self.drmaa(["gzip",  "-v",  "%s" % f], "compress")
            
    @controller.expose(help="Calculate disk usage in intermediate and data directories")
    def du(self):
        self._not_implemented()



