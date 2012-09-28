
import unittest
import tempfile
import shutil
import os
from mock import Mock
import bcbio.utils as utils

import scripts.undetermined_indexes as ui

class TestMethods(unittest.TestCase):
    
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_ui")
    
    def tearDown(self):
        shutil.rmtree(self.test_dir)
        
        
        
             
    