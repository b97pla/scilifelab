"""
Pm analysis module

Perform operations on analysis directory. 

Commands:
       ls            list contents
       runinfo       print runinfo contents
       bcstats       print information about barcode stats
       status        print status about an analysis
       hs_metrics    calculate hs metrics for samples
"""
import sys
from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.runinfo import list_runinfo

## Main analysis controller
class AnalysisController(AbstractBaseController):
    """
    Functionality for analysis management.
    """
    class Meta:
        label = 'analysis'
        description = 'Manage analysis'
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default="default")),
            (['-p', '--project'], dict(help="Project id")),
            (['-l', '--lane'], dict(help="Lane id")),
            (['-b', '--barcode_id'], dict(help="Barcode id")),
            ]

    @controller.expose(hide=True)
    def default(self):
        print __doc__

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
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        ## Get runinfo
