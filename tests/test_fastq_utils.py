"""Test suite for the fastq_utils module
"""

import tempfile
import os
import shutil
import random
import unittest
import scripts.fastq_utils as fu
import tests.generate_test_data as td

class TestFastQParser(unittest.TestCase):
    """Test the FastQParser functionality
    """
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_FastQParser_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
    def test_filter_on_header(self):
        """Filter parsing based on values in header
        """
        
        # Create some fastq records with different properties
        lanes = range(1,9)
        indexes = [td.generate_barcode() for n in xrange(4)]
        fd, fqfile = tempfile.mkstemp(suffix=".fastq.gz", dir=self.rootdir)
        os.close(fd)
        
        # For each lane and index, create a random number of records
        counts = {}
        fqw = fu.FastQWriter(fqfile)
        for l in lanes:
            counts[l] = {}
            for ix in indexes:
                no = random.randint(10,99)
                counts[l][ix] = no
                for n in xrange(no):
                    fqw.write(td.generate_fastq_record(lane=l, index=ix)) 
        fqw.close()
        
        # Try different filters and assert that the expected number of records are parsed
        
        # Filter on index
        fltr = {'index': indexes[0:-1]}
        expected = 0
        for ixc in counts.values():
            expected += sum([ixc[index] for index in fltr['index']])
        fqr = fu.FastQParser(fqfile,filter=fltr)
        for r in fqr:
            pass
        self.assertEqual(expected,fqr.rread()-1,
                         "The returned number of filtered reads based on indexes did not match expected number")

        # Filter on lane
        fltr = {'lane': range(1,5)}
        expected = 0
        for ixc in [counts[l] for l in fltr['lane']]:
            expected += sum(ixc.values())
        fqr = fu.FastQParser(fqfile,filter=fltr)
        for r in fqr:
            pass
        self.assertEqual(expected,fqr.rread()-1,
                         "The returned number of filtered reads based on lanes did not match expected number")
        
class TestFastQWriter(unittest.TestCase):
    """Test the FastQWriter functionality
    """
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_FastQWriter_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
    def test_write_fastq(self):
        """Write a fastq file
        """
        
        
        
