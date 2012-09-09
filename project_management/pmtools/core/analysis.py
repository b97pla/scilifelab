"""
Pm analysis module
"""

import sys
import os
import re
from cement.core import controller
from pmtools import AbstractBaseController
from pmtools.lib.flowcell import Flowcell
from pmtools.utils.misc import group_bcbb_files

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
            (['--pre_casava'], dict(help="Use pre-casava directory structure", action="store_true", default=False)),
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

    @controller.expose(help="Deliver data")
    def deliver(self):
        if not self._check_pargs(["flowcell", "project"]):
            return
        fc = Flowcell()
        fc.load([os.path.join(x, self.pargs.flowcell) for x in [self.config.get("archive", "root"), self.config.get("analysis", "root")]])
        if not fc:
            self.log.warn("No run information available for {}".format(self.pargs.flowcell))
            return
        indir = os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell)
        ## FIX ME: Here I'm assuming well-behaved project names
        ## collect files in project dir if casava structure
        fc.collect_files(indir)
        outdir_pfx = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "data"))
        if self.pargs.pre_casava:
            outdir = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project.replace(".", "_").lower(), "data", fc.fc_id()))
            fc_new = fc
        for sample in fc:
            key = "{}_{}".format(sample['lane'], sample['barcode_id'])
            sources = {"files":sample['files'], "results":sample['results']}
            if not self.pargs.pre_casava:
                outdir = os.path.join(outdir_pfx, sample['name'], fc.fc_id())
                fc_new = fc.subset("lane", sample['lane']).subset("name", sample['name'])
                targets = {"files": [os.path.join(outdir, os.path.basename(src)) for src in sources['files']],
                           "results": [os.path.join(outdir, os.path.basename(src)) for src in sources['results']]}
                fc_new.lane_files = dict((k, [os.path.join(outdir, os.path.basename(x)) for x in v]) for k,v in fc_new.lane_files.items())
            else:
                targets = {"files": [src.replace(indir, outdir) for src in sources['files']],
                           "results": [src.replace(indir, outdir) for src in sources['results']]}
                fc_new.lane_files = dict((k,[x.replace(indir, outdir) for x in v]) for k,v in fc_new.lane_files.items())
            fc_new.set_entry(key, 'files', targets['files'])
            fc_new.set_entry(key, 'results', targets['results'])

            if not os.path.exists(outdir):
                self.app.cmd.safe_makedir(outdir)
            ## Copy sample files - currently not doing lane files
            for src, tgt in zip(sources['files'] + sources['results'], targets['files'] + targets['results']):
                self.app.cmd.transfer_file(src, tgt)

        # for k,v in fc_new.data_dict.items():
        #     outfiles =  v['files']
        #     for f in outfiles:
        #         (pfx, bn) = (os.path.dirname(f), os.path.basename(f))
        #         # if re.match("nophix", pfx):
        #         #     continue
        #         if pfx.find("barcode"):
        #             pfx = ""
        #         if self.pargs.pre_casava:
        #             print os.path.join(outdir, pfx, bn)
        #         else:
        #             print os.path.join(outdir, v['name'], self.pargs.flowcell, f)
                
