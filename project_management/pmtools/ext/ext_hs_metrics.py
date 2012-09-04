"""
hs_metrics extension
"""
import os

import glob
from cement.core import backend, controller, handler
from pmtools import AbstractBaseController
#from pmtools.lib.runinfo import load_runinfo, subset_runinfo, runinfo_lanes, runinfo_barcodes
from pmtools.lib.flowcell import Flowcell
from pmtools.utils.misc import query_yes_no

## Auxiliary functions - move to lib or utils
def get_files(path, fc, ftype, ext=".bam", project=None, lane=None):
    """Get files from an analysis"""
    info = fc.subset("sample_prj", project)
    files = []
    for l in info.lanes():
        for bc in info.barcodes(l):
            glob_str = "{}/{}_*_*_{}-{}{}".format(path, l, bc, ftype, ext)
            glob_res = glob.glob(glob_str)
            if glob_res:
                files.append(glob_res[0])
    return files

def get_regexp_files():
    """Get files based on a regular expression in an archive folder"""
    pass

class HsMetricsController(AbstractBaseController):
    """
    Functionality for running hs_metrics
    """
    class Meta:
        label = 'hs_metrics_extension'
        """The string identifier of this handler"""
        description = 'Extension for running hs_metrics'
        stacked_on = 'analysis'
        arguments = [
            (['-r', '--region_file'], dict(help="Region definition file", default=None)),
            (['--file_type'], dict(help="File type glob", default="sort-dup")),
            ]

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["flowcell", "project", "region_file"]):
            return
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        fc = Flowcell()
        fc.load([os.path.join(x, self.pargs.flowcell) for x in [self.config.get("archive", "root"), self.config.get("analysis", "root")]])
        if not fc:
            return
        flist = get_files(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell), fc, ftype=self.pargs.file_type, project=self.pargs.project)
        if self.pargs.input_file:
            flist = [os.path.abspath(self.pargs.input_file)]
        if not query_yes_no("Going to run hs_metrics on {} files. Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        for f in flist:
            self.log.info("running CalculateHsMetrics on {}".format(f))
            ### Issue with calling java from
            ### subprocess:http://stackoverflow.com/questions/9795249/issues-with-wrapping-java-program-with-pythons-subprocess-module
            ### Actually not an issue: command line arguments have to be done the right way
            cl = ["java"] + ["-{}".format(self.pargs.java_opts)] +  ["-jar", "{}/CalculateHsMetrics.jar".format(os.getenv("PICARD_HOME"))] + ["INPUT={}".format(f)] + ["TARGET_INTERVALS={}".format(self.pargs.region_file)] + ["BAIT_INTERVALS={}".format(self.pargs.region_file)] +  ["OUTPUT={}".format(f.replace(".bam", ".hs_metrics"))] + ["VALIDATION_STRINGENCY=SILENT"]
            out = self.app.cmd.command(cl)
            if out:
                self.app._output_data["stdout"].write(out.rstrip())

def load():
    handler.register(HsMetricsController)
