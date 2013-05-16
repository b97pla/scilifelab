import os
import re
import unittest
import logbook

import subprocess 

from scilifelab.utils.misc import walk, filtered_walk, safe_makedir

filedir = os.path.abspath(__file__)
LOG = logbook.Logger(__name__)

class TestMisc(unittest.TestCase):
    def setUp(self):
        self.pattern = "^file1"

    @classmethod
    def setUpClass(self):
        dirs = ["data", "data/alignments", "data/nophix", "data/fastqc", "data/fastqc/nophix", "data/nophix/fastqc"]
        [safe_makedir(x) for x in dirs]
        [subprocess.check_call(["touch", os.path.join(x, "file1.txt")]) for x in dirs]
        [subprocess.check_call(["touch", os.path.join(x, "file2.txt")]) for x in dirs]

    def filter_fn(self, f):
        return re.search(self.pattern, f) != None

    def test_filtered_walk(self):
        """Perform a filtered walk of data dir"""
        flist = filtered_walk("data", filter_fn=self.filter_fn)
        self.assertEqual(set(flist), set(['data/file1.txt', 'data/alignments/file1.txt', 'data/nophix/file1.txt', 'data/nophix/fastqc/file1.txt', 'data/fastqc/file1.txt', 'data/fastqc/nophix/file1.txt']))

    def test_filtered_walk_include(self):
        """Perform a filtered walk of data dir, using include_dirs restriction"""
        self.pattern = "file2.txt"
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"])
        self.assertEqual(set(flist), set(['data/nophix/file2.txt', 'data/nophix/fastqc/file2.txt', 'data/fastqc/nophix/file2.txt']))

    def test_filtered_walk_exclude(self):
        """Perform a filtered walk of data dir, using exclude_dirs restriction"""
        flist = filtered_walk("data", filter_fn=self.filter_fn, exclude_dirs=["nophix"])
        self.assertEqual(set(flist), set(['data/file1.txt', 'data/alignments/file1.txt', 'data/fastqc/file1.txt']))

    def test_filtered_walk_include_exclude(self):
        """Perform a filtered walk of data dir, using include_dirs and exclude_dirs restriction"""
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"], exclude_dirs=["fastqc"])
        self.assertEqual(set(flist), set(['data/nophix/file1.txt']))

    # FIX ME: what should filtered walk return?!? If not subdirectories of fastqc then ok
    def test_filtered_walk_get_dirs(self):
        """Perform a filtered walk of data dir, getting dirs"""
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"], exclude_dirs=["fastqc"], get_dirs=True)
        self.assertEqual(set(flist), set([]))
        flist = filtered_walk("data", filter_fn=self.filter_fn, include_dirs=["nophix"], exclude_dirs=["fastqc"], get_dirs=False)
        self.assertEqual(set(flist), set(['data/nophix/file1.txt']))

