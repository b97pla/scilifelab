"""
Pm archive module

Perform operations on archive directory. 

Commands:
       ls       list contents
       runinfo  print runinfo contents
"""

import os
import yaml

from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.runinfo import runinfo_to_tab, runinfo_dump, runinfo_projects

## Main archive controller
class ArchiveController(AbstractBaseController):
    """
    Functionality for archive management.
    """
    class Meta:
        label = 'archive'
        description = 'Manage archive'
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default="default")),
            (['-p', '--project'], dict(help="Project id")),
            (['-t', '--tab'], dict(action="store_true", default=False, help="list yaml as tab file")),
            (['-P', '--list-projects'], dict(action="store_true", default=False, help="list projects of flowcell")),
            ]

    @controller.expose(hide=True)
    def default(self):
        self._not_implemented("Add help message")

    @controller.expose(help="List contents")
    def ls(self):
        self._ls("archive", "root")

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        """List runinfo for a given flowcell"""
        if self.pargs.flowcell is None or self.pargs.flowcell == "default":
            print "Please provide flowcell id"
            return
        assert self.config.get("archive", "root"), "archive directory not defined"
        f = os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml")
        self.log.info("Opening file %s" %f)
        with open(f) as fh:
            runinfo_yaml = yaml.load(fh)
        runinfo_tab = runinfo_to_tab(runinfo_yaml)
        if self.pargs.tab:
            runinfo_dump(runinfo_tab)
        elif self.pargs.list_projects:
            print "available projects for flowcell %s:\n\t" %self.pargs.flowcell + "\n\t".join(runinfo_projects(runinfo_tab))
        else:
            print runinfo_yaml

