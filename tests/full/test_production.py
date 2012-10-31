import shutil
import os
import logbook
import re
import yaml
import unittest
import drmaa

from ..classes import SciLifeTest
from classes import PmFullTest

from cement.core import handler
from scilifelab.pm.core.production import ProductionController
from scilifelab.pm.ext.ext_distributed import make_job_template_args, opt_to_dict
from scilifelab.utils.misc import filtered_walk
from scilifelab.bcbio.run import find_samples, setup_sample, remove_files, run_bcbb_command

LOG = logbook.Logger(__name__)

j_doe_00_01 = os.path.abspath(os.path.join(os.curdir, "data", "production", "J.Doe_00_01"))
j_doe_00_04 = os.path.abspath(os.path.join(os.curdir, "data", "production", "J.Doe_00_04"))
j_doe_00_05 = os.path.abspath(os.path.join(os.curdir, "data", "production", "J.Doe_00_05"))

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

ANALYSIS_TYPE = 'Align_standard_seqcap'
SAMPLE = 'P001_102_index6'
FLOWCELL = '120924_AC003CCCXX'

@unittest.skipIf(not os.getenv("MAILTO"), "not running production test: set $MAILTO environment variable to your mail address to test mailsend")
@unittest.skipIf(not os.getenv("DRMAA_LIBRARY_PATH"), "not running production test: no $DRMAA_LIBRARY_PATH")
class ProductionTest(PmFullTest):
    @classmethod
    def setUpClass(cls):
        if not os.getcwd() == filedir:
            os.chdir(filedir)
        LOG.info("Copy tree {} to {}".format(j_doe_00_01, j_doe_00_04))
        if not os.path.exists(j_doe_00_04):
            shutil.copytree(j_doe_00_01, j_doe_00_04)
        ## Set P001_102_index6 to use devel partition and require mailto environment variable for test
        pp = os.path.join(j_doe_00_04, SAMPLE, FLOWCELL, "{}-post_process.yaml".format(SAMPLE))
        with open(pp) as fh:
            config = yaml.load(fh)
        platform_args = config["distributed"]["platform_args"].split()
        platform_args[platform_args.index("-t") + 1] = "00:10:00"
        if not "--mail-user={}".format(os.getenv("MAILTO")) in platform_args:
            platform_args.extend(["--mail-user={}".format(os.getenv("MAILTO"))])
        if not "--mail-type=ALL" in platform_args:
            platform_args.extend(["--mail-type=ALL"])
        config["distributed"]["platform_args"] = " ".join(platform_args)
        with open(pp, "w") as fh:
            fh.write(yaml.safe_dump(config, default_flow_style=False, allow_unicode=True, width=1000))

    ## FIXME: since we're submitting jobs to drmaa, data will be
    ## removed before the pipeline has finished. One solution would be
    ## to run on one of the module production datasets
    # @classmethod
    # def tearDownClass(cls):
    #     LOG.info("Removing directory tree {}".format(j_doe_00_04))
    #     os.chdir(filedir)
    #     shutil.rmtree(j_doe_00_04)

    def test_production_setup(self):
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--only_setup', '--restart', '--drmaa'], extensions = ['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_production(self):
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--amplicon', '--restart'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_platform_args(self):
        """Test the platform arguments for a run"""
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--amplicon', '--restart', '--sample', SAMPLE, '--drmaa'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_change_platform_args(self):
        """Test that passing --time actually changes platform
        arguments. These arguments should have precedence over
        whatever is written in the config file."""
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--amplicon', '--restart', '--sample', SAMPLE, '--drmaa', '--time', '00:01:00', '-n'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_casava_transfer(self):
        """Test transfer of casava data from production to project"""
        self.app = self.make_app(argv = ['production', 'transfer', 'J.Doe_00_03', '--debug', '--force', '--quiet'], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)
        j_doe_00_03 = os.path.abspath(os.path.join(os.curdir, "data", "projects", "j_doe_00_03"))
        pattern = ".fastq(.gz)?$"
        def fastq_filter(f):
            return re.search(pattern, f) != None
        fastq_files = filtered_walk(j_doe_00_03, fastq_filter)
        self.assertEqual(len(fastq_files), 2)

    def test_touch_finished(self):
        """Test touching finished files"""
        SAMPLE2 = "P001_101_index3"
        finished = os.path.join(os.curdir, "data", "production", "J.Doe_00_01", SAMPLE, "FINISHED_AND_DELIVERED")
        finished2 = os.path.join(os.curdir, "data", "production", "J.Doe_00_01", SAMPLE2, "FINISHED_AND_DELIVERED")
        if os.path.exists(finished):
            os.unlink(finished)
        if os.path.exists(finished2):
            os.unlink(finished2)
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_01', '--debug', '--force', '--sample', SAMPLE], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        self.assertTrue(os.path.exists(finished))
        samplefile = os.path.join(os.curdir, "data", "production", "J.Doe_00_01", "finished_sample.txt")
        with open(samplefile, "w") as fh:
            fh.write(SAMPLE + "\n")
            fh.write(SAMPLE2 + "\n")
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_01', '--debug', '--force', '--sample', samplefile], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        self.assertTrue(os.path.exists(finished2))
        ## Make sure rsync fails
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_01', '--debug', '--force', '--sample', samplefile], extensions=[])
        handler.register(ProductionController)
        try:
            self.app.setup()
            self.app.config.set("runqc", "root", self.app.config.get("runqc", "root").replace("production", "projects"))
            with self.app.log.log_setup.applicationbound():
                self.app.run()
                self.app.render(self.app._output_data)
        finally:
            self.app.close()

class UtilsTest(SciLifeTest):
    @classmethod
    def setUpClass(cls):
        if not os.getcwd() == filedir:
            os.chdir(filedir)
        LOG.info("Copy tree {} to {}".format(j_doe_00_01, j_doe_00_05))
        if not os.path.exists(j_doe_00_05):
            shutil.copytree(j_doe_00_01, j_doe_00_05)

    @classmethod
    def tearDownClass(cls):
        LOG.info("Removing directory tree {}".format(j_doe_00_05))
        os.chdir(filedir)
        shutil.rmtree(j_doe_00_05)

    def test_find_samples(self):
        """Test finding samples"""
        flist = find_samples(j_doe_00_05)
        self.assertEqual(len(flist), 3)
        flist = find_samples(j_doe_00_05, **{'only_failed':True})
        self.assertEqual(len(flist), 0)

    def test_setup_samples(self):
        """Test setting up samples, changing genome to rn4"""
        flist = find_samples(j_doe_00_05)
        for f in flist:
            setup_sample(f, **{'analysis_type':'Align_standard_seqcap', 'genome_build':'rn4', 'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'num_cores':8, 'distributed':False})
        for f in flist:
            with open(f, "r") as fh:
                config = yaml.load(fh)
            self.assertEqual(config["details"][0]["multiplex"][0]["genome_build"], "rn4")

            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            self.assertIn("--no-google-report", cl)
            self.assertIn("--only-run", cl)
            with open(f.replace("-bcbb-config.yaml", "-post_process.yaml")) as fh:
                config = yaml.load(fh)
            self.assertEqual(config["custom_algorithms"][ANALYSIS_TYPE]["hybrid_bait"], 'rat_baits.interval_list')
            self.assertEqual(config["custom_algorithms"][ANALYSIS_TYPE]["hybrid_target"], 'rat_targets.interval_list')
            self.assertEqual(config["algorithm"]["num_cores"], 8)
                

        for f in flist:
            setup_sample(f, **{'analysis_type':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':True, 'google_report':True, 'analysis_type':'Align_standard_seqcap'
                               , 'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':False})
            with open(f, "r") as fh:
                config = yaml.load(fh)
            self.assertEqual(config["details"][0]["multiplex"][0]["genome_build"], "rn4")
            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            self.assertNotIn("--no-google-report", cl)
            self.assertNotIn("--only-run", cl)
            with open(f.replace("-bcbb-config.yaml", "-post_process.yaml")) as fh:
                config = yaml.load(fh)
            self.assertEqual(config["algorithm"]["mark_duplicates"], False)
            self.assertEqual(config["custom_algorithms"][ANALYSIS_TYPE]["mark_duplicates"], False)

    def test_remove_files(self):
        """Test removing files"""
        keep_files = ["-post_process.yaml$", "-post_process.yaml.bak$", "-bcbb-config.yaml$", "-bcbb-config.yaml.bak$",  "-bcbb-command.txt$", "-bcbb-command.txt.bak$", "_[0-9]+.fastq$", "_[0-9]+.fastq.gz$", "^[0-9][0-9]_.*.txt$"]
        pattern = "|".join(keep_files)
        def remove_filter_fn(f):
            return re.search(pattern, f) == None
        flist = find_samples(j_doe_00_05)
        for f in flist:
            workdir = os.path.dirname(f)
            remove_files = filtered_walk(workdir, remove_filter_fn)
            self.assertNotIn("01_analysis_start.txt", [os.path.basename(x) for x in remove_files])

    def test_remove_dirs(self):
        """Test removing directories before rerunning pipeline"""
        keep_files = ["-post_process.yaml$", "-post_process.yaml.bak$", "-bcbb-config.yaml$", "-bcbb-config.yaml.bak$",  "-bcbb-command.txt$", "-bcbb-command.txt.bak$", "_[0-9]+.fastq$", "_[0-9]+.fastq.gz$"]
        pattern = "|".join(keep_files)
        def remove_filter_fn(f):
            return re.search(pattern, f) == None
        flist = find_samples(j_doe_00_05)
        for f in flist:
            workdir = os.path.dirname(f)
            remove_dirs = filtered_walk(workdir, remove_filter_fn, get_dirs=True)
            self.assertIn("fastqc", [os.path.basename(x) for x in remove_dirs])
        
    def test_bcbb_command(self):
        """Test output from command, changing analysis to amplicon and
        setting targets and baits"""
        flist = find_samples(j_doe_00_05)
        for f in flist:
            setup_sample(f, **{'analysis_type':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':False, 'google_report':False, 'analysis_type':'Align_standard_seqcap',
                               'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':False})
            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            (cl, platform_args) = run_bcbb_command(f)
            self.assertIn("automated_initial_analysis.py",cl)
            setup_sample(f, **{'analysis_type':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':False, 'google_report':False, 'analysis_type':'Align_standard_seqcap',
                               'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':True})
            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            (cl, platform_args) = run_bcbb_command(f)
            self.assertIn("distributed_nextgen_pipeline.py",cl)
    
    @unittest.skipIf(not os.getenv("DRMAA_LIBRARY_PATH"), "not running UtilsTest.test_platform: no $DRMAA_LIBRARY_PATH")
    def test_platform_args(self):
        """Test making platform args and changing them on the fly. """
        pp = os.path.join(j_doe_00_05, SAMPLE, FLOWCELL, "{}-post_process.yaml".format(SAMPLE))
        with open(pp) as fh:
            config = yaml.load(fh)
        platform_args = config["distributed"]["platform_args"].split()
        self.assertIn("core", platform_args)
        pargs = opt_to_dict(platform_args)
        self.assertEqual("P001_102_index6-bcbb.log", pargs['-o'])
        kw = {'time':'00:01:00', 'jobname':'test', 'partition':'devel'}
        pargs = make_job_template_args(pargs, **kw)
        self.assertEqual("devel", pargs['partition'])
        nativeSpec = "-t {time} -p {partition} -A {account}".format(**pargs)
        self.assertEqual("00:01:00", nativeSpec[3:11])
                    

