"""
Pm compress module

Functions for compressing files in archive, project and analysis folders.

"""
import os
import sys

from cement.core import controller, output, backend
from cement.utils.shell import *
from pmtools import AbstractBaseController, SubSubController

Log = backend.minimal_logger(__name__)

## Main compress controller
class CompressController(AbstractBaseController):
    """
    Functionality for project management.
    """
    class Meta:
        label = 'compress'
        stacked_on = None
        description = 'compress utilities'
        arguments = [
            (['id'], dict(help="project/flowcell id", default=None, action="store")),
            (['-f', '--fastq'], dict(help="compress fastq files. Implicitly sets pattern ", default=False, action="store_true")),
            (['-d', '--decompress'], dict(help="decompress", default=False, action="store_true")),
            (['-s', '--sbatch'], dict(help="run as sbatch", default=False, action="store_true")),
            (['--pattern'], dict(help="pattern to search for", default=None, action="store")),
            ]

    @controller.expose(hide=True)
    def default(self):
        print __doc__

    @controller.expose(help="perform compress operation on a project in project folder")
    def project(self):
        if not os.path.exists(os.path.join(self.config.get("projects", "root"), self.pargs.id)):
            self.app.log.warn("No such folder: %s" % os.path.join(self.config.get("projects", "root"), self.pargs.id))
            sys.exit()

    @controller.expose(help="perform compress operation on a flowcell in analysis folder")
    def analysis(self):
        self._not_implemented()


    
