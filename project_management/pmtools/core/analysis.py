"""
Pm analysis module
"""

import sys
import os
from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.flowcell import Flowcell
from pmtools.utils.misc import group_bcbb_files

## Main analysis controller
class AnalysisController(AbstractBaseController):
    """
    Functionality for analysis management.
    """
    class Meta:
        label = 'analysis'
        description = 'Manage analysis'
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default=None)),
            (['-p', '--project'], dict(help="Project id")),
            (['-l', '--lane'], dict(help="Lane id")),
            (['-b', '--barcode_id'], dict(help="Barcode id")),
            ]

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="List contents")
    def ls(self):
        self._ls("analysis", "root")

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        self._not_implemented()

    @controller.expose(help="List bcstats")
    def bcstats(self):
        self._not_implemented()

    @controller.expose(help="List status of a run")
    def status(self):
        self._not_implemented()

    @controller.expose(help="Deliver data")
    def deliver(self):
        if not self._check_pargs(["flowcell", "project"]):
            return
        fc = Flowcell()
        fc.load([os.path.join(x, self.pargs.flowcell) for x in [self.config.get("archive", "root"), self.config.get("analysis", "root")]])
        if not fc:
            return
        flist = fc.get_files(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell), ftype=self.pargs.file_type, project=self.pargs.project)
        ## FIX ME: Here I'm assuming well-behaved project names
        outdir = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "data", self.pargs.flowcell))
        print dir(fc)
        print fc.dict()
        # fdictlist = [group_bcbb_files(x) for x in flist]
        # print fdictlist
        
