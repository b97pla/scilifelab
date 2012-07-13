"""
Pm analysis module

Perform operations on analysis directory. 

Commands:
       ls       list contents
       runinfo  print runinfo contents
       bcstats  print information about barcode stats
       status   print status about an analysis
"""

from cement.core import controller
from pmtools import AbstractBaseController

## Main analysis controller
class AnalysisController(AbstractBaseController):
    """
    Functionality for analysis management.
    """
    class Meta:
        label = 'analysis'
        description = 'Manage analysis'

    @controller.expose(hide=True)
    def default(self):
        print __doc__
        sys.exit()

    @controller.expose(help="List contents")
    def ls(self):
        self._not_implemented()

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        self._not_implemented()

    @controller.expose(help="List bcstats")
    def bcstats(self):
        self._not_implemented()

    @controller.expose(help="List status of a run")
    def status(self):
        self._not_implemented()
