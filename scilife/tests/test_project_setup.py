"""
Test project setup

Once project fastq files have been delivered with sample_delivery.py they have 
to be relinked to comply with bcbio.
"""
import os
import sys
import subprocess
import unittest
import shutil

class ProjectSetupTest(unittest.TestCase):
    """Test project setup"""

    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_01")

    def test_project_setup(self):
        """Test project setup"""
        cl = ["setup_project_files.py",
              os.path.join(self.proj_dir, "project_run_info.yaml"),
              "20000101A_hiseq2000",
              "--project_dir=%s" %(self.proj_dir)]
        print cl
        subprocess.check_call(cl)
