"""
Pm clean module

Utilities for removing files.

FIXME: - rename to rm?
       - more help
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
            (['id'], dict(help="Flowcell id/project id", action="store", default=None)),
            (['-l', '--lane'], dict(help="Lane id", default=None, action="store")),
            (['-p', '--project'], dict(help="Project id", action="store")),
            (['-b', '--bam_files'], dict(help="Work on bam files", default=False, action="store_true")),
            (['--pattern'], dict(help="pattern to clean", default=None, action="store")),
            ]

    def _setup(self, app_obj):
        # shortcuts
        super(AbstractBaseController, self)._setup(app_obj)
    
    @controller.expose(help="Clean up folder", hide=True)
    def default(self):
        """Clean up folder"""
        print self.__doc__

    @controller.expose(help="perform clean operation on a project in project folder")
    def project(self):
        ## Set the pattern to use
        pattern = None
        if self.pargs.bam_files:
            pattern = ".bam$"
        else:
            pattern = self.pargs.pattern
        if not pattern:
            self.log.info("No pattern defined")
            return

        def keep(f):
            return re.search(pattern, f) != None
        flist = []
        for root, dirs, files in os.walk(os.path.join(self.config.get("projects", "root"), "projects")):
            flist = flist + filter(keep, files)
                  
    @controller.expose(help="perform clean operation on a flowcell in analysis folder")
    def analysis(self):
        self._not_implemented()
        print "FIXME: implement Mayas code here"

        
