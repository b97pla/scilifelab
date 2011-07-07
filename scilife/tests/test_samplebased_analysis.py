"""This directory is setup with configurations to run the main functional test.

It exercises a samplebased analysis pipeline on a smaller subset of data, as implemented at SciLife.
"""
import os
import sys
import subprocess
import unittest
import shutil
import contextlib
import glob

try:
    import bcbio
except:
    raise ImportError("Module bcbio required to run sample based analysis")

@contextlib.contextmanager
def make_workdir():
    dirname = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "intermediate")
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)

class SampleBasedAnalysisTest(unittest.TestCase):
    """Setup a sample based scilife analysis
    """
    def setUp(self):
        self.proj_dir = os.path.join(os.pardir, "projects", "j_doe_00_01")

    def _data_delivery(self):
        """Deliver data from bcbio"""
        bcbio_tests_dir = os.path.join(bcbio.__path__[0], os.pardir, "tests")
        bcbio_data_dir = os.path.join(bcbio_tests_dir, "data", "110106_FC70BUKAAXX")

        if not os.path.exists(bcbio_data_dir):
            raise IOError("bcbio example data has not been installed: you need to run nosetests in your bcbio test directory " + bcbio_tests_dir)

        data_delivery_dirs = {'20000101A_hiseq2000':os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "data", "20000101A_hiseq2000"),
                              '20000101B_hiseq2000':os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_01", "data", "20000101B_hiseq2000")}
        run_samples = {'20000101A_hiseq2000' : {'1_110106_FC70BUKAAXX_1_fastq.txt' : '1_000101_FC70BUKAAXX_SAMPLE1_index1_1.fastq', 
                                                '1_110106_FC70BUKAAXX_2_fastq.txt' : '1_000101_FC70BUKAAXX_SAMPLE1_index1_2.fastq', 
                                                '2_110106_FC70BUKAAXX_1_fastq.txt' : '2_000101_FC70BUKAAXX_SAMPLE1_index2_1.fastq',
                                                '2_110106_FC70BUKAAXX_2_fastq.txt' : '2_000101_FC70BUKAAXX_SAMPLE1_index2_2.fastq',
                                                '3_110106_FC70BUKAAXX_fastq.txt' : '3_000101_FC70BUKAAXX_SAMPLE2_index3.fastq'},
                       '20000101B_hiseq2000' : {'1_110106_FC70BUKAAXX_1_fastq.txt' : '1_000101_FC70BUKAAXX_SAMPLE1_index1_1.fastq', 
                                                '1_110106_FC70BUKAAXX_2_fastq.txt' : '1_000101_FC70BUKAAXX_SAMPLE1_index1_2.fastq', 
                                                '2_110106_FC70BUKAAXX_1_fastq.txt' : '2_000101_FC70BUKAAXX_SAMPLE2_index2_1.fastq',
                                                '2_110106_FC70BUKAAXX_2_fastq.txt' : '2_000101_FC70BUKAAXX_SAMPLE2_index2_2.fastq'}
                       }
        for name, dirname in data_delivery_dirs.items():
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            samples = run_samples[name]
            for src, dest in samples.items():
                infile = os.path.join(bcbio_data_dir, src)
                outfile = os.path.join(dirname, dest)
                if not os.path.exists(outfile):
                    print "Coyping " + infile + " to " +  outfile
                    shutil.copyfile(infile, outfile)
                else:
                    print "Data already delivered for " + outfile + " in " + dirname

    def _install_exome_target_files(self):
        """Download and install target test files"""
        pass

    def test_run_samplebased_pipeline(self):
        self._data_delivery()
        script_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "scripts")

        with make_workdir():
            cl = [os.path.join(script_dir, "exome_pipeline.py"),
                  os.path.join(self.proj_dir, "post_process.yaml"),
                  os.path.join(self.proj_dir, "data", "20000101A_hiseq")]
            subprocess.check_call(cl)
