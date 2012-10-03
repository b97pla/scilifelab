"""Test suite for the fastq_utils module
"""

import tempfile
import os
import shutil
import random
import unittest
import copy
import scilifelab.utils.fastq_utils as fu
import tests.generate_test_data as td
import scilifelab.illumina.hiseq as hi

class TestFastQParser(unittest.TestCase):
    """Test the FastQParser functionality
    """
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_FastQParser_")
        
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
        
        self.example_fq = fqfile
        self.example_counts = counts
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
    def test_filter_on_header(self):
        """Filter parsing based on values in header
        """
        
        # Try different filters and assert that the expected number of records are parsed
        
        # Filter on index
        fltr = {'index': self.example_counts.values()[0].keys()[0:-1]}
        expected = 0
        for ixc in self.example_counts.values():
            expected += sum([ixc[index] for index in fltr['index']])
        fqr = fu.FastQParser(self.example_fq,filter=fltr)
        for r in fqr:
            pass
        self.assertEqual(expected,fqr.rread()-1,
                         "The returned number of filtered reads based on indexes did not match expected number")

        # Filter on lane
        fltr = {'lane': range(1,5)}
        expected = 0
        for ixc in [self.example_counts[l] for l in fltr['lane']]:
            expected += sum(ixc.values())
        fqr = fu.FastQParser(self.example_fq,filter=fltr)
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
        

class TestFastQUtils(unittest.TestCase):
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_fastq_utils_")
        
        # Generate some test fastq_files
        indexes = {td.generate_barcode(): [],
                   td.generate_barcode(): [],
                   td.generate_barcode(): [],
                   td.generate_barcode(): []}
        args = {'instrument': td.generate_instrument(),
                'run_number': random.randint(101,9999),
                'fcid': td.generate_fc_barcode(),
                'lane': 1,
                'pair': True}
        fd, f1 = tempfile.mkstemp(suffix="_R1.fastq.gz", dir=self.rootdir)
        os.close(fd) 
        f2 = f1.replace("_R1","_R2")
        f1h = fu.FastQWriter(f1)
        f2h = fu.FastQWriter(f2)
        for n in range(1000):
            args['index'] = random.choice(indexes.keys()) 
            record = td.generate_fastq_record(**args)
            indexes[args['index']].append(record[0])
            f1h.write(record[0:4])
            f2h.write(record[4:])
        f1h.close()
        f2h.close()
        
        # Create a samplesheet to use for demultiplexing using all but the last index
        samplesheet = f1.replace("_R1.fastq.gz",".csv")
        sdata = []
        for n, index in enumerate(indexes.keys()[0:-1]):
             sdata.append([args['fcid'],
                           str(args['lane']),
                           "Sample_{}".format(str(n)),
                           "unknown",
                           index,
                           "DemuxTest",
                           "0",
                           "",
                           "",
                           "DemuxTestProject"])
             
        # Add an entry that will not match
        s = copy.copy(sdata[-1])
        s[1] = str(args['lane']+1)
        sdata.append(s)
        
        self.samplesheet = td._write_samplesheet(sdata,samplesheet)
        self.indexes = indexes
        self.fastq_1 = f1
        self.fastq_2 = f2
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
             
    def test_demultiplex_fastq(self):
        """Demultiplex a test fastq file
        """
        
        # Demultiplex sample files based on samplesheet
        outfiles = fu.demultiplex_fastq(self.rootdir,self.samplesheet,self.fastq_1,self.fastq_2)
        
        # Assert that no files were returned for empty lanes
        n = 0
        for lane in outfiles.keys():
            n += sum([len(ix) for ix in outfiles[lane].values() if lane != "1"])
        self.assertEqual(0,
                         n,
                         "Demultiplexing should not return results for empty lane")
        
        # Assert that the expected number of output files were returned
        outfiles = outfiles["1"]
        sdata = hi.HiSeqRun.parse_samplesheet(self.samplesheet)
        self.assertEqual(len([s for s in sdata if s["Lane"] == "1"]),
                         len(outfiles.keys()),
                         "Demultiplexing did not return the expected number of fastq files")
        
        for index in outfiles.keys():
            # Assert that the out_files was written to the correct folder
            self.assertEqual([self.rootdir,self.rootdir],
                             [os.path.dirname(o) for o in outfiles[index]],
                             "The demultiplexed output was not written to the correct folder")
            
            # Parse the outfile and verify the output
            headers = []
            f1h = fu.FastQParser(outfiles[index][0])
            f2h = fu.FastQParser(outfiles[index][1])
            for r1 in f1h:
                r2 = f2h.next()
                r1s = r1[0].strip().split()
                r2s = r2[0].strip().split()
                self.assertListEqual([r1s[0],r1s[1][1:]],
                                     [r2s[0],r2s[1][1:]],
                                     "Header strings from paired fastq files don't match")
                headers.append(r1[0])
            
            # Assert that the number of sequences matches the expected and that the headers match
            self.assertEqual(len(headers),len(self.indexes[index]),
                             "The number of demultiplexed reads in file does not match expected")
            self.assertListEqual(sorted(headers),sorted(self.indexes[index]),
                                 "The parsed headers from demultiplexed fastq file do not match the expected")
            
   
        
