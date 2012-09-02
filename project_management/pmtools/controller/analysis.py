"""Pm analysis module"""

usage = """Pm analysis module

Perform operations on analysis directory. 

Commands:
       ls            list contents
       runinfo       print runinfo contents
       bcstats       print information about barcode stats
       status        print status about an analysis
       hs_metrics    calculate hs metrics for samples
"""
import sys
import os
from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.runinfo import get_runinfo, subset_runinfo

## Auxiliary functions
def get_files(runinfo_tab, type="fastq", project=None, lane=None):
    """Get files from an analysis"""
    i = subset_runinfo(runinfo_tab, "sample_prj", project)
    return i

def get_regexp_files():
    """Get files based on a regular expression in an archive folder"""
    pass

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
        print usage

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

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["flowcell", "project"]):
            return
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        ## Get runinfo
        if os.path.exists(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml")):
            runinfo_tab = get_runinfo(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml"))
        elif os.path.exists(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell, "run_info.yaml")):
            runinfo_tab = get_runinfo(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell, "run_info.yaml"))
        else:
            self.log.warn("No run information available")
            return
        print get_files(runinfo_tab, project=self.pargs.project)
