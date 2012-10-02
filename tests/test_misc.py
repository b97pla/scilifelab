import os
import re
import unittest

import subprocess 

from scilifelab.utils.misc import *

filedir = os.path.abspath(__file__)

class TestMisc(unittest.TestCase):
    def setUp(self):
        dirs = ["data", "data/alignments", "data/nophix", "data/fastqc", "data/fastqc/nophix", "data/nophix/fastqc"]
        [safe_makedir(x) for x in dirs]
        [subprocess.check_call(["touch", os.path.join(x, "file1.txt")]) for x in dirs]
        [subprocess.check_call(["touch", os.path.join(x, "file2.txt")]) for x in dirs]
        self.pattern = "^file1"
        
    def filter_fn(self, f):
        return re.search(self.pattern, f) != None

    def test_1_filtered_walk(self):
        """Perform a filtered walk of data dir"""
        flist = filtered_walk("data", filter_fn=self.filter_fn)
        self.assertEqual(flist, ['data/file1.txt', 'data/alignments/file1.txt', 'data/nophix/file1.txt', 'data/nophix/fastqc/file1.txt', 'data/fastqc/file1.txt', 'data/fastqc/nophix/file1.txt'])


    def test_2_filtered_walk_include(self):
        """Perform a filtered walk of data dir, using include_dirs restriction"""
        self.pattern = "file2.txt"
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"])
        self.assertEqual(flist, ['data/nophix/file2.txt', 'data/nophix/fastqc/file2.txt', 'data/fastqc/nophix/file2.txt'])


    def test_3_filtered_walk_exclude(self):
        """Perform a filtered walk of data dir, using exclude_dirs restriction"""
        flist = filtered_walk("data", filter_fn=self.filter_fn, exclude_dirs=["nophix"])
        self.assertEqual(flist, ['data/file1.txt', 'data/alignments/file1.txt', 'data/fastqc/file1.txt'])


    def test_4_filtered_walk_include_exclude(self):
        """Perform a filtered walk of data dir, using include_dirs and exclude_dirs restriction"""
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"], exclude_dirs=["fastqc"])
        self.assertEqual(flist, ['data/nophix/file1.txt'])

