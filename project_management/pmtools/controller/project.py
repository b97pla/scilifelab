"""
Pm Project module

Provide functionality for project management.

Commands:
       ls          list contents
       init        initialize a project folder
       add         add boilerplate code
       compress    compress files
       clean       remove files 
"""

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

    @controller.expose(hide=True)
    def default(self):
        __doc__

    @controller.expose(help="List project folder")
    def ls(self):
        self._not_implemented()

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
        self._not_implemented()



