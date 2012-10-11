"""hs_metrics extension"""
import os
import re

import glob
from cement.core import backend, controller, handler
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.bcbio.flowcell import Flowcell
from scilifelab.utils.misc import query_yes_no, filtered_walk

class HsMetricsController(AbstractBaseController):
    """
    Functionality for running hs_metrics
    """
    class Meta:
        label = 'hs_metrics_extension'
        """The string identifier of this handler"""
        description = 'Extension for running hs_metrics'
        stacked_on = 'production'
        arguments = [
            (['--region_file'], dict(help="Region definition file", default=None)),
            (['--bait_file'], dict(help="Region bait definition file", default=None)),
            ## FIX ME: This should be called bcbb_file_type and be loaded via an extension
            (['--hs_file_type'], dict(help="File type glob", default="sort-dup")),
            ]

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["project", "region_file"]):
            return
        if not self.pargs.bait_file:
            self.pargs.bait_file = self.pargs.region_file
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        pattern = "{}.bam$".format(self.pargs.hs_file_type)
        def filter_fn(f):
            return re.search(pattern, f) != None
        ### FIX ME: this isn't caught by _process_args
        path =  self.pargs.flowcell if self.pargs.flowcell else self.pargs.project
        flist = filtered_walk(os.path.join(self.config.get("production", "root"), path), filter_fn=filter_fn, exclude_dirs=['nophix', 'alignments', 'fastqc', 'fastq_screen'])
        if self.pargs.input_file:
            flist = [os.path.abspath(self.pargs.input_file)]
        if not query_yes_no("Going to run hs_metrics on {} files. Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        for f in flist:
            self.log.info("running CalculateHsMetrics on {}".format(f))
            ### Issue with calling java from
            ### subprocess:http://stackoverflow.com/questions/9795249/issues-with-wrapping-java-program-with-pythons-subprocess-module
            ### Actually not an issue: command line arguments have to be done the right way
            cl = ["java"] + ["-{}".format(self.pargs.java_opts)] +  ["-jar", "{}/CalculateHsMetrics.jar".format(os.getenv("PICARD_HOME"))] + ["INPUT={}".format(f)] + ["TARGET_INTERVALS={}".format(os.path.abspath(self.pargs.region_file))] + ["BAIT_INTERVALS={}".format(os.path.abspath(self.pargs.bait_file))] +  ["OUTPUT={}".format(f.replace(".bam", ".hs_metrics"))] + ["VALIDATION_STRINGENCY=SILENT"]
            out = self.app.cmd.command(cl)
            if out:
                self.app._output_data["stdout"].write(out.rstrip())

def load():
    handler.register(HsMetricsController)
