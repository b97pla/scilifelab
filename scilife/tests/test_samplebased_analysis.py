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
from string import Template

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
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")

    def _data_delivery(self):
        """Deliver data from bcbio"""
        bcbio_tests_dir = os.path.join(bcbio.__path__[0], os.pardir, "tests")
        bcbio_barcode_dirs = {'1' : os.path.join(bcbio_tests_dir, "test_automated_output", "1_110106_FC70BUKAAXX_barcode"),
                              '2' : os.path.join(bcbio_tests_dir, "test_automated_output", "2_110106_FC70BUKAAXX_barcode"),
                              '3' : os.path.join(bcbio_tests_dir, "test_automated_output", "3_110106_FC70BUKAAXX_barcode")
                              }

        if not os.path.exists(bcbio_barcode_dirs['1']):
            raise IOError("bcbio example data has not been installed: you need to run nosetests in your bcbio test directory " + bcbio_tests_dir)

        data_delivery_dirs = {'20000101A_hiseq2000':os.path.join(self.file_dir, "projects", "j_doe_00_01", "data", "20000101A_hiseq2000"),
                              '20000101B_hiseq2000':os.path.join(self.file_dir, "projects", "j_doe_00_01", "data", "20000101B_hiseq2000")}
        run_samples = {'20000101A_hiseq2000' : {'1_110106_FC70BUKAAXX_2_1_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE1_index10_1.fastq',
                                                '1_110106_FC70BUKAAXX_2_2_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE1_index10_2.fastq',
                                                '1_110106_FC70BUKAAXX_3_1_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE1_index2_1.fastq',
                                                '1_110106_FC70BUKAAXX_3_2_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE1_index2_2.fastq',
                                                '1_110106_FC70BUKAAXX_4_1_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE2_index3_1.fastq',
                                                '1_110106_FC70BUKAAXX_4_2_fastq.txt' : '1_110106_FC70BUKAAXX_SAMPLE2_index3_2.fastq',
                                                '2_110106_FC70BUKAAXX_6_1_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE1_index1_1.fastq',
                                                '2_110106_FC70BUKAAXX_6_2_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE1_index1_2.fastq',
                                                '2_110106_FC70BUKAAXX_7_1_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE2_index4_1.fastq',
                                                '2_110106_FC70BUKAAXX_7_2_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE2_index4_2.fastq'
                                                },
                       '20000101B_hiseq2000' : {'2_110106_FC70BUKAAXX_8_1_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE1_index3_1.fastq',
                                                '2_110106_FC70BUKAAXX_8_2_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE1_index3_2.fastq',
                                                '2_110106_FC70BUKAAXX_9_1_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE2_index4_1.fastq',
                                                '2_110106_FC70BUKAAXX_9_2_fastq.txt' : '2_110106_FC70BUKAAXX_SAMPLE2_index4_2.fastq'
                                                }
                       }
        for name, dirname in data_delivery_dirs.items():
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            samples = run_samples[name]
            for src, dest in samples.items():
                barcode = src.split("_")[0]
                infile = os.path.join(bcbio_barcode_dirs[barcode], src)
                outfile = os.path.join(dirname, dest)
                if not os.path.exists(outfile):
                    print "Coyping " + infile + " to " +  outfile
                    shutil.copyfile(infile, outfile)
                else:
                    print "Data already delivered as " + os.path.basename(outfile) + " in " + dirname

        loc_files = ['bowtie_indices.loc', 'bwa_index.loc', 'sam_fa_indices.loc']
        tooldir = os.path.join(self.file_dir, "config", "tool-data")

        d = {'genomedir':os.path.join(bcbio_tests_dir, "data", "genomes")}
        if not os.path.exists(tooldir):
            os.makedirs(tooldir)
        for lf in loc_files:
            outfile = os.path.join(tooldir, lf)
            if not os.path.exists(outfile):
                with open(os.path.join(self.file_dir, "templates", "tool-data", lf)) as in_handle:
                    tmpl = Template(in_handle.read())
                with open(outfile, "w") as out_handle:
                    out_handle.write(tmpl.safe_substitute(d))
                          

    def _install_exome_target_files(self):
        """Download and install target test files"""
        pass

    def _install_project_config_files(self): 
        """Install the project config files, inserting the correct paths to galaxy configuration files"""
        proj_conf = os.path.join(self.file_dir, "templates", "proj_conf.yaml")
        d = {'galaxy_config' : os.path.join(self.file_dir, "config", "universe_wsgi.ini"),
             'log_dir' : os.path.join(self.file_dir, "log")}
        with open(proj_conf) as in_handle:
            tmpl = Template(in_handle.read())
        print "Installing project configuration in " + self.proj_dir
        with open (os.path.join(self.proj_dir, "proj_conf.yaml"), "w") as out_handle:
            out_handle.write(tmpl.safe_substitute(d))

    def test_run_samplebased_pipeline(self):
        self._data_delivery()
        self._install_project_config_files()
        script_dir = os.path.join(self.file_dir, os.pardir, os.pardir, "scripts")

        with make_workdir():
            cl = [os.path.join(script_dir, "exome_pipeline.py"),
                  os.path.join(self.proj_dir, "proj_conf.yaml"),
                  os.path.join(self.proj_dir, "data", "20000101A_hiseq2000")]
            subprocess.check_call(cl)
