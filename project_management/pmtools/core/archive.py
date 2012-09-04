"""Pm archive module"""

import os
import yaml

from cement.core import controller
from pmtools import AbstractBaseController
#from pmtools.lib.runinfo import get_runinfo, runinfo_projects, dump_runinfo
from pmtools.lib.flowcell import *

## Main archive controller
class ArchiveController(AbstractBaseController):
    """
    Functionality for archive management.

    This is the base controller for archive management.
    """
    class Meta:
        """Controller meta-data settings"""

        label = 'archive'
        description = 'Manage archive'
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default="default")),
            (['-p', '--project'], dict(help="Project id")),
            (['--as_yaml'], dict(action="store_true", default=False, help="list runinfo as yaml file")),
            (['-P', '--list-projects'], dict(action="store_true", default=False, help="list projects of flowcell")),
            ]

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="List contents")
    def ls(self):
        return self._ls("archive", "root")

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        """List runinfo for a given flowcell"""
        if self.pargs.flowcell is None or self.pargs.flowcell == "default":
            self.app._output_data["stderr"].write("Please provide flowcell id")
            return
        assert self.config.get("archive", "root"), "archive directory not defined"
        fc = Flowcell(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml"))
        self.log.info("Opening file {}".format(fc.filename))
        if self.pargs.list_projects:
            self.app._output_data['stdout'].write("available projects for flowcell {}:\n\t".format(self.pargs.flowcell) + "\n\t".join(fc.projects()))
            return
        if self.pargs.as_yaml:
            self.app._output_data['stdout'].write(str(fc.as_yaml()))
        else:            
            self.app._output_data['stdout'].write(str(fc))

