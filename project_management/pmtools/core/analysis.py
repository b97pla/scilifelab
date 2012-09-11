"""
Pm analysis module
"""

import sys
import os
import re
from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.flowcell import Flowcell

## Main analysis controller
class AnalysisController(AbstractBaseController):
    """
    Functionality for analysis management.
    """
    class Meta:
        label = 'analysis'
        description = 'Manage analysis'
        arguments = [
            (['flowcell'], dict(help="Flowcell id", nargs="?", default=None)),
            (['-p', '--project'], dict(help="Project id")),
            (['-l', '--lane'], dict(help="Lane id")),
            (['-b', '--barcode_id'], dict(help="Barcode id")),
            (['--from_pre_casava'], dict(help="Use pre-casava directory structure for gathering information", action="store_true", default=False)),
            (['--to_pre_casava'], dict(help="Use pre-casava directory structure for delivery", action="store_true", default=False)),
            ]

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="List contents")
    def ls(self):
        self._ls(self.app.config.get("analysis", "root"))

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        self._not_implemented()

    @controller.expose(help="List bcstats")
    def bcstats(self):
        self._not_implemented()

    @controller.expose(help="List status of a run")
    def status(self):
        self._not_implemented()


    def _from_casava_structure(self):
        """Get information from casava structure"""
        if not self._check_pargs(["project"]):
            return
        self._not_implemented("Reading of casava structure not yet tested nor supported")

    def _to_casava_structure(self, fc):
        outdir_pfx = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "data"))
        for sample in fc:
            key = "{}_{}".format(sample['lane'], sample['barcode_id'])
            sources = {"files":sample['files'], "results":sample['results']}
            outdir = os.path.join(outdir_pfx, sample['name'], fc.fc_id())
            dirs = {"data":os.path.abspath(os.path.join(outdir_pfx, sample['name'], fc.fc_id())),
                    "intermediate":os.path.abspath(os.path.join(outdir_pfx, sample['name'], fc.fc_id()))}
            self._make_output_dirs(dirs)
            fc_new = fc.subset("lane", sample['lane']).subset("name", sample['name'])
            targets = {"files": [os.path.join(dirs["data"], os.path.basename(src)) for src in sources['files']],
                       "results": [os.path.join(dirs["intermediate"], os.path.basename(src)) for src in sources['results']]}
            fc_new.lane_files = dict((k, [os.path.join(outdir, os.path.basename(x)) for x in v]) for k,v in fc_new.lane_files.items())
            fc_new.set_entry(key, 'files', targets['files'])
            fc_new.set_entry(key, 'results', targets['results'])
            ## Copy sample files - currently not doing lane files
            self._transfer_files(sources, targets)

    def _to_pre_casava_structure(self, fc):
        dirs = {"data":os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "data", fc.fc_id())),
                "intermediate":os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "intermediate", fc.fc_id()))}
        self._make_output_dirs(dirs)
        fc_new = fc
        for sample in fc:
            key = "{}_{}".format(sample['lane'], sample['barcode_id'])
            sources = {"files":sample['files'], "results":sample['results']}
            targets = {"files": [src.replace(fc.path, dirs["data"]) for src in sources['files']],
                       "results": [src.replace(fc.path, dirs["intermediate"]) for src in sources['results']]}
            fc_new.set_entry(key, 'files', targets['files'])
            fc_new.set_entry(key, 'results', targets['results'])
            ## FIX ME: lane file gathering
            ## fc_new.lane_files = dict((k,[x.replace(indir, outdir) for x in v]) for k,v in fc_new.lane_files.items())
            ## Copy sample files - currently not doing lane files
            self._transfer_files(sources, targets)
        with open(os.path.join(dirs["data"], "project_run_info.yaml"), "w") as yaml_out:
            yaml_out.write(fc_new.as_yaml())

    def _from_pre_casava_structure(self):
        if not self._check_pargs(["project", "flowcell"]):
            return
        fc = Flowcell()
        fc.load([os.path.join(x, self.pargs.flowcell) for x in [self.config.get("archive", "root"), self.config.get("analysis", "root")]])
        indir = os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell)
        if not fc:
            self.log.warn("No run information available for {}".format(self.pargs.flowcell))
            return
        fc_new = fc.subset("sample_prj", self.pargs.project)
        fc_new.collect_files(indir)        
        return fc_new

    def _make_output_dirs(self, dirs):
        if not os.path.exists(dirs["data"]):
            self.app.cmd.safe_makedir(dirs["data"])
        if not os.path.exists(dirs["intermediate"]):
            self.app.cmd.safe_makedir(dirs["intermediate"])

    def _transfer_files(self, sources, targets):
        for src, tgt in zip(sources['files'] + sources['results'], targets['files'] + targets['results']):
            if not os.path.exists(os.path.dirname(tgt)):
                self.app.cmd.safe_makedir(os.path.dirname(tgt))
            self.app.cmd.transfer_file(src, tgt)


    @controller.expose(help="Deliver data")
    def deliver(self):
        if not self.pargs.from_pre_casava and self.pargs.to_pre_casava:
            self.app.log.warn("not delivering from casava input to pre_casava output")
            return
        ## Collect files depending on input structure
        if self.pargs.from_pre_casava:
            fc = self._from_pre_casava_structure()
        else:
            fc = self._from_casava_structure()
        if not fc:
            return
        ## Organize output file names depending on output structure
        if self.pargs.to_pre_casava:
            self._to_pre_casava_structure(fc)
        else:
            self._to_casava_structure(fc)

        ## Fix lane_files for pre_casava output
        ## if self.pargs.pre_casava:

