import shutil
import os
import logbook
import re
import yaml
import unittest
import pandas as pd
import numpy as np

try:
    import drmaa
except:
    pass

from ..classes import SciLifeTest
from classes import PmFullTest

from cement.core import handler
from scilifelab.pm.core.production import ProductionController
from scilifelab.utils.misc import filtered_walk, opt_to_dict
from scilifelab.bcbio.run import find_samples, setup_sample, run_bcbb_command, setup_merged_samples, sample_table, get_vcf_files, validate_sample_directories

LOG = logbook.Logger(__name__)

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

j_doe_00_01 = os.path.abspath(os.path.join(filedir, "data", "production", "J.Doe_00_01"))
j_doe_00_04 = os.path.abspath(os.path.join(filedir, "data", "production", "J.Doe_00_04"))
j_doe_00_05 = os.path.abspath(os.path.join(filedir, "data", "production", "J.Doe_00_05"))

ANALYSIS_TYPE = 'Align_standard_seqcap'
GALAXY_CONFIG = os.path.abspath(os.path.join(filedir, "data", "config"))
SAMPLES =  ['P001_101_index3', 'P001_102_index6']
FLOWCELL = '120924_AC003CCCXX'

FINISHED = {
    'J.Doe_00_01': {'P001_101_index3': os.path.join(filedir, "data", "production", "J.Doe_00_01", SAMPLES[0], "FINISHED_AND_DELIVERED"),
                    'P001_102_index6': os.path.join(filedir, "data", "production", "J.Doe_00_01", SAMPLES[1], "FINISHED_AND_DELIVERED")},
    'J.Doe_00_04': {'P001_101_index3': os.path.join(filedir, "data", "production", "J.Doe_00_04", SAMPLES[0], "FINISHED_AND_DELIVERED"),
                    'P001_102_index6': os.path.join(filedir, "data", "production", "J.Doe_00_04", SAMPLES[1], "FINISHED_AND_DELIVERED")}
    }

REMOVED = {
    'J.Doe_00_01': {'P001_101_index3': os.path.join(filedir, "data", "production", "J.Doe_00_01", SAMPLES[0], "FINISHED_AND_REMOVED"),
                    'P001_102_index6': os.path.join(filedir, "data", "production", "J.Doe_00_01", SAMPLES[1], "FINISHED_AND_REMOVED")},
    'J.Doe_00_04': {'P001_101_index3': os.path.join(filedir, "data", "production", "J.Doe_00_04", SAMPLES[0], "FINISHED_AND_REMOVED"),
                    'P001_102_index6': os.path.join(filedir, "data", "production", "J.Doe_00_04", SAMPLES[1], "FINISHED_AND_REMOVED")}
    }



class ProductionConsoleTest(PmFullTest):
    """Class for testing functions without drmaa"""

    def setUp(self):
        pass

    def test_compression_suite(self):
        """Test various combinations of compression, decompression, cleaning"""
        self.app = self.make_app(argv = ['production', 'decompress', 'J.Doe_00_01', '--debug', '--force', '--fastq',  '-n'])
        handler.register(ProductionController)
        self._run_app()
        l1 = self.app._output_data["stderr"].getvalue()
        self.app = self.make_app(argv = ['production', 'decompress', 'J.Doe_00_01', '-f', FLOWCELL, '--debug', '--force', '--fastq',  '-n'])
        handler.register(ProductionController)
        self._run_app()
        l2 = self.app._output_data["stderr"].getvalue()
        self.assertTrue(len(l1) > len(l2))
        os.chdir(filedir)

    def test_run(self):
        """Test various combinations of compression, decompression, cleaning"""
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_01', '--debug', '--force', '--fastq',  '-n'])
        handler.register(ProductionController)
        self._run_app()
        l1 = self.app._output_data["stderr"].getvalue()
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_01', '-f', FLOWCELL, '--debug', '--force', '--fastq',  '-n'])
        handler.register(ProductionController)
        self._run_app()
        l2 = self.app._output_data["stderr"].getvalue()
        self.assertTrue(len(l1) > len(l2))
        os.chdir(filedir)
        


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
        pp = os.path.join(j_doe_00_04, SAMPLES[1], FLOWCELL, "{}-post_process.yaml".format(SAMPLES[1]))
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
        for k in FINISHED.keys():
            for v in FINISHED[k].values():
                if os.path.exists(v):
                    os.unlink(v)
            for v in REMOVED[k].values():
                if os.path.exists(v):
                    os.unlink(v)

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
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--amplicon', '--restart', '--sample', SAMPLES[1], '--drmaa'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_change_platform_args(self):
        """Test that passing --time actually changes platform
        arguments. These arguments should have precedence over
        whatever is written in the config file."""
        self.app = self.make_app(argv = ['production', 'run', 'J.Doe_00_04', '--debug', '--force', '--amplicon', '--restart', '--sample', SAMPLES[1], '--drmaa', '--time', '00:01:00', '-n'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)

    def test_casava_transfer(self):
        """Test transfer of casava data from production to project"""
        self.app = self.make_app(argv = ['production', 'transfer', 'J.Doe_00_03', '--debug', '--force', '--quiet'], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        os.chdir(filedir)
        j_doe_00_03 = os.path.abspath(os.path.join(filedir, "data", "projects", "j_doe_00_03"))
        pattern = ".fastq(.gz)?$"
        def fastq_filter(f):
            return re.search(pattern, f) != None
        fastq_files = filtered_walk(j_doe_00_03, fastq_filter)
        self.assertEqual(len(fastq_files), 2)

    def test_touch_finished(self):
        """Test touching finished files"""
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_01', '--debug', '--force', '--sample', SAMPLES[0]], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        self.assertTrue(os.path.exists(FINISHED['J.Doe_00_01'][SAMPLES[0]]))
        samplefile = os.path.join(filedir, "data", "production", "J.Doe_00_01", "finished_sample.txt")
        with open(samplefile, "w") as fh:
            fh.write(SAMPLES[0] + "\n")
            fh.write(SAMPLES[1] + "\n")
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_01', '--debug', '--force', '--sample', samplefile], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        self.assertTrue(os.path.exists(FINISHED['J.Doe_00_01'][SAMPLES[1]]))
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

    def test_remove_finished(self):
        self.app = self.make_app(argv = ['production', 'touch-finished', 'J.Doe_00_04', '--debug', '--force', '--sample', SAMPLES[1]], extensions=[])
        handler.register(ProductionController)
        self._run_app()
        self.assertTrue(os.path.exists(FINISHED['J.Doe_00_04'][SAMPLES[1]]))
        ## Remove file, dry
        self.app = self.make_app(argv = ['production', 'remove-finished', 'J.Doe_00_04', '--debug', '--force', '-n'], extensions=[])
        handler.register(ProductionController)
        self._run_app()



class UtilsTest(SciLifeTest):
    @classmethod
    def setUpClass(cls):
        if not os.getcwd() == filedir:
            os.chdir(filedir)
        LOG.info("Copy tree {} to {}".format(j_doe_00_01, j_doe_00_05))
        if not os.path.exists(j_doe_00_05):
            shutil.copytree(j_doe_00_01, j_doe_00_05)
            with open(os.path.join(j_doe_00_05, "samples.txt"), "w") as fh:
                fh.write("\n\nP001_101_index3\nP001_104_index3")
            with open(os.path.join(j_doe_00_05, "samples2.txt"), "w") as fh:
                fh.write("\n\nP001_101_index3-bcbb-config.yaml")

    @classmethod
    def tearDownClass(cls):
        LOG.info("Removing directory tree {}".format(j_doe_00_05))
        os.chdir(filedir)
        shutil.rmtree(j_doe_00_05)

    def test_find_samples(self):
        """Test finding samples"""
        flist = find_samples(j_doe_00_05)
        self.assertIn(len(flist), [3,4])
        flist = find_samples(j_doe_00_05, **{'only_failed':True})
        self.assertIn(len(flist), [0,1])

    def test_find_samples_from_file(self):
        """Find samples defined in file with empty lines and erroneous names"""
        with open(os.path.join(j_doe_00_05, "P001_101_index3-bcbb-config.yaml"), "w") as fh:
            fh.write("\n")
        flist = find_samples(j_doe_00_05, sample=os.path.join(j_doe_00_05, "samples.txt"))
        validate_sample_directories(flist, j_doe_00_05)
        self.assertEqual(len(flist),2)
        os.unlink(os.path.join(j_doe_00_05, "P001_101_index3-bcbb-config.yaml"))

    def test_find_samples_from_file_with_yaml(self):
        """Find samples defined in file with empty lines and a bcbb-config.yaml file lying directly under root directory"""
        flist = find_samples(j_doe_00_05, sample=os.path.join(j_doe_00_05, "samples2.txt"))
        args = [flist, j_doe_00_05]
        self.assertRaises(Exception, validate_sample_directories, *args)

    def test_setup_merged_samples(self):
        """Test setting up merged samples"""
        flist = find_samples(j_doe_00_05)
        setup_merged_samples(flist, **{'dry_run':False})

    def test_setup_samples(self):
        """Test setting up samples, changing genome to rn4"""
        flist = find_samples(j_doe_00_05)
        for f in flist:
            setup_sample(f, **{'analysis':'Align_standard_seqcap', 'genome_build':'rn4', 'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'num_cores':8, 'distributed':False})
        for f in flist:
            with open(f, "r") as fh:
                config = yaml.load(fh)
            if config["details"][0].get("multiplex", None):
                self.assertEqual(config["details"][0]["multiplex"][0]["genome_build"], "rn4")
            else:
                self.assertEqual(config["details"][0]["genome_build"], "rn4")

            with open(f.replace("-bcbb-config.yaml", "-post_process.yaml")) as fh:
                config = yaml.load(fh)
            self.assertEqual(config["custom_algorithms"][ANALYSIS_TYPE]["hybrid_bait"], 'rat_baits.interval_list')
            self.assertEqual(config["custom_algorithms"][ANALYSIS_TYPE]["hybrid_target"], 'rat_targets.interval_list')
            self.assertEqual(config["algorithm"]["num_cores"], 8)
                
        for f in flist:
            setup_sample(f, **{'analysis':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':True, 'google_report':True,
                               'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':False})
            with open(f, "r") as fh:
                config = yaml.load(fh)
            if config["details"][0].get("multiplex", None):
                self.assertEqual(config["details"][0]["multiplex"][0]["genome_build"], "rn4")
            else:
                self.assertEqual(config["details"][0]["genome_build"], "rn4")
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
            setup_sample(f, **{'analysis':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':False, 'google_report':False,
                               'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':False})
            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            (cl, platform_args) = run_bcbb_command(f)
            self.assertIn("automated_initial_analysis.py",cl)
            setup_sample(f, **{'analysis':ANALYSIS_TYPE, 'genome_build':'rn4', 'dry_run':False,
                               'no_only_run':False, 'google_report':False, 
                               'dry_run':False, 'baits':'rat_baits.interval_list', 'targets':'rat_targets.interval_list', 'amplicon':True, 'num_cores':8, 'distributed':True})
            with open(f.replace("-bcbb-config.yaml", "-bcbb-command.txt")) as fh:
                cl = fh.read().split()
            (cl, platform_args) = run_bcbb_command(f)
            self.assertIn("distributed_nextgen_pipeline.py",cl)

    def test_global_post_process(self):
        """Test that when using a "global" post_process, jobname,
        output, error and output directory are updated.
        """
        flist = find_samples(j_doe_00_05)
        pp = os.path.join(j_doe_00_01, SAMPLES[1], FLOWCELL, "{}-post_process.yaml".format(SAMPLES[1]))
        with open(pp) as fh:
            postprocess = yaml.load(fh)
        for f in flist:
            (cl, platform_args) = run_bcbb_command(f, pp)
            self.assertIn("--error", platform_args)
            self.assertEqual(platform_args[platform_args.index("--error") + 1], f.replace("-bcbb-config.yaml", "-bcbb.err"))
    
    @unittest.skipIf(not os.getenv("DRMAA_LIBRARY_PATH"), "not running UtilsTest.test_platform: no $DRMAA_LIBRARY_PATH")
    def test_platform_args(self):
        """Test making platform args and changing them on the fly. """
        from scilifelab.pm.ext.ext_distributed import make_job_template_args
        pp = os.path.join(j_doe_00_05, SAMPLES[1], FLOWCELL, "{}-post_process.yaml".format(SAMPLES[1]))
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
                    

    def test_sample_table(self):
        """Test making a sample table"""
        flist = find_samples(j_doe_00_01)
        samples = sample_table(flist)
        grouped = samples.groupby("sample")
        self.assertEqual(len(grouped.groups["P001_101_index3"]), 2)
        self.assertEqual(len(grouped.groups["P001_102_index6"]), 1)

    def test_summarize_variants(self):
        """Test summarizing variants"""
        flist = find_samples(j_doe_00_01)
        vcf_d = get_vcf_files(flist)
            
