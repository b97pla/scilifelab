"""
Pm ls module

Functions for listing archive, analysis and project directories
"""
import re
import sys
import os
import yaml
import csv

from cement.core import controller, output, backend
from cement.utils.shell import *
from pmtools import AbstractBaseController, SubSubController

Log = backend.minimal_logger(__name__)

## Main ls controller
class LsController(AbstractBaseController):
    """
    Functionality for listing archive, analysis, project folders
    """
    class Meta:
        label = 'ls'
        stacked_on = None
        description = 'List archive, analysis, project folders'
        interface = controller.IController
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default="default")),
            (['-a', '--archive'], dict(help="List archive")),
            (['-p', '--project'], dict(help="Project id")),
            ]

    @controller.expose(hide=True)
    def default(self):
        print "FIXME: show help"

    def _setup(self, app_obj):
        # shortcuts
        super(AbstractBaseController, self)._setup(app_obj)
        # Compile ignore regexps
        self.reignore = re.compile(self.config.get("config", "ignore").replace("\n", "|"))

    def _filtered_ls(self, out):
        """Filter output"""
        def ignore(line):
            return self.reignore.match(line) == None
        return filter(ignore, out)
    
    @controller.expose(help="List project folder")
    def proj(self):
        assert self.config.get("projects", "root"), "no projects root directory"
        (out, err, code) = exec_cmd(["ls",  self.config.get("projects", "root")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)

    @controller.expose(help="List finished projects folder")
    def finished_proj(self):
        assert self.config.get("projects","root"), "no projects root directory"
        (out, err, code) = exec_cmd(["ls",  os.path.join(self.config.get("projects", "root"), "finished_projects")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)

    @controller.expose(help="List archive folder")
    def archive(self):
        """List contents of archive folder"""
        assert self.config.get("config", "archive"), "no config archive directory"
        (out, err, code) = exec_cmd(["ls",  self.config.get("config", "archive")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)

    @controller.expose(help="List analysis folder")
    def analysis(self):
        """List contents of analysis folder"""
        (out, err, code) = exec_cmd(["ls",  self.app.config.get("config", "analysis")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)

## Runinfo controller
class RunInfoController(SubSubController):
    class Meta:
        label = 'runinfo-ctrl'
        stacked_on= 'ls'
        description='Runinfo controller'
        interface = controller.IController
        arguments = [
            ##(["flowcell"], dict(type=str)),
            (['-t', '--tab'], dict(action="store_true", default=False, help="list yaml as tab file")),
            (['-P', '--projects'], dict(action="store_true", default=False, help="list projects of yaml file")),
            ]

    def _yaml_to_tab(self, runinfo_yaml):
        """Convert yaml to tab"""
        out = []
        lanelabels = ["lane", "description", "flowcell_id", "genome_build", "analysis"]
        mplabels = ["barcode_type", "barcode_id", "sample_prj", "name", "sequence"]
        header = lanelabels + mplabels
        out.append(header)
        for info in runinfo_yaml:
            laneinfo = [info.get(x, None) for x in lanelabels]
            for mp in info.get("multiplex", None):
                mpinfo = [mp.get(x, None) for x in mplabels]
                line = laneinfo + mpinfo
                out.append(line)
        return out

    def _column(self, matrix, label):
        i = matrix[0].index(label)
        return [row[i] for row in matrix[1:]]

    @controller.expose(hide=True)
    def default(self):
        pass

    @controller.expose(help="List runinfo file for a given flowcell")
    def runinfo(self):
        """List runinfo for a given flowcell"""
        if self.pargs.flowcell is None:
            print "Please provide flowcell id"
            return
        assert self.config.get("config", "archive"), "archive directory not defined"
        f = os.path.join(self.config.get("config", "archive"), self.pargs.flowcell, "run_info.yaml")
        self.log.info("Opening file %s" %f)
        with open(f) as fh:
            runinfo_yaml = yaml.load(fh)
        runinfo_tab = self._yaml_to_tab(runinfo_yaml)
        if self.pargs.tab:
            w=csv.writer(sys.stdout, delimiter="\t")
            w.writerows(runinfo_tab)
        elif self.pargs.projects:
            print "available projects for flowcell %s:\n\t" %self.pargs.flowcell + "\n\t".join(list(set(self._column(runinfo_tab, "sample_prj"))))
        else:
            print runinfo_yaml
