"""
Pm clean module
"""
import re
import sys
import os
import yaml

from cement.core import controller, output, backend
from cement.utils.shell import *
from pmtools import AbstractBaseController, SubSubController

Log = backend.minimal_logger(__name__)

## Main production controller
class CleanController(AbstractBaseController):
    """Provides basic functions for cleaning up analysis/project folders"""
    class Meta:
        label = 'clean'
        stacked_on = None
        description = 'Functionality for cleaning up analysis/project folders'
        interface = controller.IController
        arguments = [
            (['fc'], dict(help="flowcell", action="store", default=None, nargs="?")),
            (['-f', '--flowcell'], dict(help="Flowcell id", action="store")),
            (['-p', '--project'], dict(help="Project id", action="store")),
            (['-b', '--bam_files'], dict(help="Work on bam files", default=False, action="store_true")),
            ]

    def _setup(self, app_obj):
        # shortcuts
        super(AbstractBaseController, self)._setup(app_obj)
    
    @controller.expose(help="Clean up folder", hide=True)
    def default(self):
        """Clean up folder"""
        if self.pargs.fc is None:
            print self.__doc__
        else:
            print self.pargs

    
