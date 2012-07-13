"""
Pm archive module

Perform operations on archive directory. 

Commands:
       ls       list contents
       runinfo  print runinfo contents
"""

from cement.core import controller
from pmtools import AbstractBaseController

## Main archive controller
class ArchiveController(AbstractBaseController):
    """
    Functionality for archive management.
    """
    class Meta:
        label = 'archive'
        description = 'Manage archive'

    @controller.expose(hide=True)
    def default(self):
        __doc__

    @controller.expose(help="List contents")
    def ls(self):
        self._not_implemented()
