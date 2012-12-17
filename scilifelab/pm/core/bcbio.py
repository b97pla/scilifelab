"""pm bcbio extension"""
import os
import re
import yaml
import shutil

from cement.core import backend, controller, handler, hook
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.misc import query_yes_no, filtered_walk
from scilifelab.bcbio.run import find_samples, setup_sample, remove_files, run_bcbb_command, setup_merged_samples, get_vcf_files, validate_sample_directories
from scilifelab.report.qc import compile_qc 
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

BCBIO_EXCLUDE_DIRS = ['realign-split', 'variants-split', 'tmp', 'tx', 'fastqc', 'fastq_screen', 'alignments', 'nophix']

class BcbioRunController(AbstractBaseController):
    class Meta:
        label = 'bcbio'
        description = 'Wrapper for bcbio analyses'

    def _setup(self, app):
        group = app.args.add_argument_group('Bcbio argument group', 'Options for bcbio')
        group.add_argument('--post_process', help="post process file. Setting this will override sample-specific post process files.", action="store", default=None, nargs="?", type=str)
        group.add_argument('--genome_build', help="genome build ", action="store", default=None, type=str)
        group.add_argument('--only_failed', help="only run on failed samples ", action="store_true", default=False)
        group.add_argument('--amplicon', help="amplicon-based analyses (e.g. HaloPlex), which means mark_duplicates is set to false", action="store_true", default=False)
        group.add_argument('--targets', help="sequence capture target file", action="store", default=None)
        group.add_argument('--baits', help="sequence capture baits file", action="store", default=None)
        group.add_argument('--distributed', help="run distributed, changing 'num_cores' in  post_process to 'messaging': calls automated_initial_analysis.py", action="store_true", default=False)
        group.add_argument('--num_cores', help="num_cores value; default 8", action="store", default=8, type=int)
        group.add_argument('--only_setup', help="only perform setup", action="store_true", default=False)
        group.add_argument('--restart', help="restart analysis", action="store_true", default=False)
        group.add_argument('--analysis', help="set analysis type in bcbb config file", action="store", default=None, type=str)
        group.add_argument('--snpEff', help="set snpEff program in post process", action="store", default=None, type=str)
        group.add_argument('--merged', help="do merged sample analysis for samples with multiple sample runs", action="store_true", default=False)
        group.add_argument('--new_config', help="make new config file", action="store_true", default=False)
        group.add_argument('--hs_file_type', help="File type glob", default="sort-dup")
        super(BcbioRunController, self)._setup(app)

    def _sample_status(self, x):
        """Find the status of a sample.

        Look for output files: currently only look for project-summary.csv"""
        if os.path.exists(os.path.join(os.path.dirname(x), "project-summary.csv")):
            return "PASS"
        else:
            return "FAIL"

    @controller.expose(help="Run bcbb pipeline")
    def run(self):
        if not self._check_pargs(["project"]):
            return
        if self.pargs.post_process:
            self.pargs.post_process = os.path.abspath(self.pargs.post_process)
        basedir = os.path.abspath(os.path.join(self.app.controller._meta.root_path, self.app.controller._meta.path_id))
        flist = find_samples(basedir, **vars(self.pargs))
        # Add filtering on flowcell if necessary
        self._meta.pattern = ".*"
        flist = [x for x in flist if self._filter_fn(x)]
        if self.pargs.merged:
            ##  Setup merged samples and append to flist if new list longer
            flist = setup_merged_samples(flist, **vars(self.pargs))
        if not len(flist) > 0:
            self.log.info("No sample configuration files found")
            return
        if len(flist) > 0 and not query_yes_no("Going to start {} jobs... Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        # Make absolutely sure analysis directory is a *subdirectory* of the working directory
        validate_sample_directories(flist, basedir)
        orig_dir = os.path.abspath(os.getcwd())

        for run_info in flist:
            os.chdir(os.path.abspath(os.path.dirname(run_info)))
            setup_sample(run_info, **vars(self.pargs))
            os.chdir(orig_dir)
        if self.pargs.only_setup:
            return
        if self.pargs.only_failed:
            status = {x:self._sample_status(x) for x in flist}
            flist = [x for x in flist if self._sample_status(x)=="FAIL"]
        ## Here process files again, removing if requested, and running the pipeline
        for run_info in flist:
            self.app.log.info("Running analysis defined by config file {}".format(run_info))
            os.chdir(os.path.abspath(os.path.dirname(run_info)))
            if self.app.cmd.monitor(work_dir=os.path.dirname(run_info)):
                self.app.log.warn("Not running job")
                continue
            if self.pargs.restart:
                self.app.log.info("Removing old analysis files in {}".format(os.path.dirname(run_info)))
                remove_files(run_info, **vars(self.pargs))
            (cl, platform_args) = run_bcbb_command(run_info, **vars(self.pargs))
            self.app.cmd.command(cl, **{'platform_args':platform_args, 'saveJobId':True, 'workingDirectory':os.path.dirname(run_info)})
            os.chdir(orig_dir)

    @controller.expose(help="Compile qc metrics based on result files")
    def compile_qc(self):
        """Compile qc metrics for samples without statusdb information."""
        if not self._check_pargs(["project"]):
            return
        kw = {'exclude_dirs': BCBIO_EXCLUDE_DIRS}
        kw.update(**vars(self.pargs))
        out_data = compile_qc(path=os.path.abspath(os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id)), **kw)
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())

    @controller.expose(help="Calculate hs metrics for samples")
    def hs_metrics(self):
        if not self._check_pargs(["project", "targets"]):
            return
        if not self.pargs.baits:
            self.pargs.baits = self.pargs.targets
        self.log.info("hs_metrics: This is a temporary solution for calculating hs metrics for samples using picard tools")
        pattern = "{}.bam$".format(self.pargs.hs_file_type)
        def filter_fn(f):
            return re.search(pattern, f) != None
        ### FIX ME: this isn't caught by _process_args
        flist = []
        path =  self.pargs.flowcell if self.pargs.flowcell else self.pargs.project
        if self.pargs.sample:
            if os.path.exists(self.pargs.sample):
                with open(self.pargs.sample) as fh:
                    flist = [x.rstrip() for x in fh.readlines()]
            else:
                pattern = "{}{}".format(sample, pattern)
        if not flist:
            flist = filtered_walk(os.path.join(self.config.get(self.app.controller._meta.label, "root"), path), filter_fn=filter_fn, exclude_dirs=['nophix', 'alignments', 'fastqc', 'fastq_screen'])
        if not query_yes_no("Going to run hs_metrics on {} files. Are you sure you want to continue?".format(len(flist)), force=self.pargs.force):
            return
        for f in flist:
            self.log.info("running CalculateHsMetrics on {}".format(f))
            ### Issue with calling java from
            ### subprocess:http://stackoverflow.com/questions/9795249/issues-with-wrapping-java-program-with-pythons-subprocess-module
            ### Actually not an issue: command line arguments have to be done the right way
            cl = ["java"] + ["-{}".format(self.pargs.java_opts)] +  ["-jar", "{}/CalculateHsMetrics.jar".format(os.getenv("PICARD_HOME"))] + ["INPUT={}".format(f)] + ["TARGET_INTERVALS={}".format(os.path.abspath(self.pargs.targets))] + ["BAIT_INTERVALS={}".format(os.path.abspath(self.pargs.baits))] +  ["OUTPUT={}".format(f.replace(".bam", ".hs_metrics"))] + ["VALIDATION_STRINGENCY=SILENT"]
            out = self.app.cmd.command(cl)
            if out:
                self.app._output_data["stdout"].write(out.rstrip())

    @controller.expose(help="Perform basic variant summary")
    def vcf_summary(self):
        if not self._check_pargs(["project"]):
            return
        flist = find_samples(os.path.abspath(os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id)), **vars(self.pargs))
        vcf_d = get_vcf_files(flist)
        ## Traverse files, copy to result directory, run bgzip and tabix, and merge vcfs to one file
        outdir = os.path.join(os.path.abspath(os.path.join(self.app.controller._meta.project_root, self.app.controller._meta.path_id, "intermediate", "results", "vcf")))
        if not os.path.exists(outdir):
            self.app.cmd.safe_makedir(outdir)
        for k, v in vcf_d.iteritems():
            print v
            if v.endswith(".gz"):
                tgt = os.path.join(outdir, os.path.basename(v).replace("TOTAL", "TOTAL_{}".format(k)))
                v = v.replace(".gz", "")
                tgt = tgt.replace(".gz", "")
            else:
                ## bgzip
                LOG.info("Running bgzip on {}".format(v))
                cl = ["bgzip", v]
                self.app.cmd.command(cl)
            ##if not os.path.exists("{}.gz.tbi"):
            ## tabix
            LOG.info("Running tabix on {}.gz".format(v))
            cl = ["tabix", "-f", "-p", "vcf", "{}.gz".format(v)]
            self.app.cmd.command(cl)
            self.app.cmd.link("{}.gz".format(v), "{}.gz".format(tgt))
            self.app.cmd.link("{}.gz.tbi".format(v), "{}.gz.tbi".format(tgt))
        ## Make all-variants file
        # all_variants = os.path.join(outdir, "all-variants.vcf")
        # cl = vcf-merge `ls -1 *.vcf.gz | tr '\n', ' '` > all-variants.vcf
