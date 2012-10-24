import os
import sys
import subprocess
import unittest
import collections

from Bio import SeqIO

from scilifelab.utils.misc import safe_makedir

class SciLifeTest(unittest.TestCase):
    """Test class for small unit tests that only test functions and do
    not rely on the existence of pipeline-generated data.
    """
    def setUp(self):
        pass

class SciLifeFullTest(unittest.TestCase):
    """Test class for tests that require that test data has been
    completely run through the pipeline. The class will setup test
    fastq files and run them through the pipeline once.
    """

    def setUp(self):
        pass

