"""
hs_metrics extension
"""
import os
import glob
from cement.core import backend, controller, handler
from pmtools import AbstractBaseController
from pmtools.lib.runinfo import load_runinfo, subset_runinfo, runinfo_lanes, runinfo_barcodes
from pmtools.utils.misc import query_yes_no

## Auxiliary functions
def get_files(path, runinfo_tab, ext="sort-dup.bam", project=None, lane=None):
    """Get files from an analysis"""
    info = subset_runinfo(runinfo_tab, "sample_prj", project)
    files = []
    lanes = runinfo_lanes(info)
    for l in lanes:
        barcodes = runinfo_barcodes(info, l)
        for bc in barcodes:
            glob_str = "{}/{}_*_*_{}-{}".format(path, l, bc, ext)
            files.append(glob.glob(glob_str)[0])
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
            ]

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["flowcell", "project", "region_file"]):
            return
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        runinfo_tab = load_runinfo(self.config, self.pargs.flowcell)
        if not runinfo_tab:
            return
        flist = get_files(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell), runinfo_tab, project=self.pargs.project)
        if not query_yes_no("Going to run hs_metrics on {} files. Are you sure you want to continue?".format(len(flist))):
            sys.exit()
        for f in flist:
            self.log.info("running CalculateHsMetrics on {}".format(f))
            out = self.app.cmd.command(["java -jar $PICARD_HOME/CalculateHsMetrics.jar -I {} -TI {} -BI {} -O {}".format(f, self.pargs.region_file, self.pargs.region_file, f.replace(".bam", ".hs_metrics"))])
            if out:
                self.app._output_data["stdout"].write(out.rstrip())

def load():
    handler.register(HsMetricsController)
