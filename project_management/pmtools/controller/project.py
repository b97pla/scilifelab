"""
Pm Project module
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

    @controller.expose()
    def init(self):
        print "in init"


