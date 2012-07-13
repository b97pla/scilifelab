"""
Pm archive module

Perform operations on archive directory. 

Commands:
       ls       list contents
       runinfo  print runinfo contents
"""

import os

from cement.core import controller
from cement.utils.shell import *
from pmtools import AbstractBaseController


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
        print __doc__


    @controller.expose(help="List contents")
    def ls(self):
        assert self.config.get("archive", "root"), "no config archive directory"
        (out, err, code) = exec_cmd(["ls",  self.config.get("archive", "root")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)

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
        elif self.pargs.projects:
            print "available projects for flowcell %s:\n\t" %self.pargs.flowcell + "\n\t".join(list(set(self._column(runinfo_tab, "sample_prj"))))
        else:
            print runinfo_yaml

