"""
hs_metrics extension
"""
from cement.core import backend, controller, handler
from pmtools import AbstractBaseController

class HsMetricsController(AbstractBaseController):
    """
    Functionality for running hs_metrics
    """
    class Meta:
        label = 'hs_metrics_extension'
        """The string identifier of this handler"""
        description = 'Extension for running hs_metrics'
        stacked_on = 'analysis'

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["flowcell", "project"]):
            return
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        ## Get runinfo
        if os.path.exists(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml")):
            runinfo_tab = get_runinfo(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml"))
        elif os.path.exists(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell, "run_info.yaml")):
            runinfo_tab = get_runinfo(os.path.join(self.config.get("analysis", "root"), self.pargs.flowcell, "run_info.yaml"))
        else:
            self.log.warn("No run information available")
            return
        print get_files(runinfo_tab, project=self.pargs.project)

def load():
    handler.register(HsMetricsController)
