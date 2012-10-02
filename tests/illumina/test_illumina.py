
import unittest
import shutil
import tempfile
import os
import random
import string
import bcbio.utils as utils
import tests.generate_test_data as td
from scripts.fastq_utils import (FastQParser, FastQWriter, demultiplex_fastq)
import scilifelab.illumina as illumina
from scilifelab.illumina import IlluminaRun

class TestIllumina(unittest.TestCase):
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_illumina_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
    
    def test_demultiplex_fastq(self):
        """Demultiplex a test fastq file
        """
        
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
        f1h = FastQWriter(f1)
        f2h = FastQWriter(f2)
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
        samplesheet = td._write_samplesheet(sdata,samplesheet)
        
        # Demultiplex sample files based on samplesheet
        outfiles = demultiplex_fastq(self.rootdir,samplesheet,f1,f2)
        outfiles = outfiles["1"]
        
        # Assert that the expected number of output files were returned
        self.assertEqual(len(sdata),len(outfiles.keys()),
                         "Demultiplexing did not return the expected number of fastq files")
        
        for index in outfiles.keys():
            # Assert that the out_files was written to the correct folder
            self.assertEqual([self.rootdir,self.rootdir],
                             [os.path.dirname(o) for o in outfiles[index]],
                             "The demultiplexed output was not written to the correct folder")
            
            # Parse the outfile and verify the output
            headers = []
            f1h = FastQParser(outfiles[index][0])
            f2h = FastQParser(outfiles[index][1])
            for r1 in f1h:
                r2 = f2h.next()
                r1s = r1[0].strip().split()
                r2s = r2[0].strip().split()
                self.assertListEqual([r1s[0],r1s[1][1:]],
                                     [r2s[0],r2s[1][1:]],
                                     "Header strings from paired fastq files don't match")
                headers.append(r1[0])
            
            # Assert that the number of sequences matches the expected and that the headers match
            self.assertEqual(len(headers),len(indexes[index]),
                             "The number of demultiplexed reads in file does not match expected ({} vs {})".format(str(len(headers)),str(len(indexes[index]))))
            self.assertListEqual(sorted(headers),sorted(indexes[index]),
                                 "The parsed headers from demultiplexed fastq file do not match the expected")
            
            
            
        
        
class TestIlluminaRun(unittest.TestCase):
    
    def setUp(self): 
        self.rootdir = tempfile.mkdtemp(prefix="test_illumina_run_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
    
    def test__get_flowcell(self):
        """Get flowcell from analysis directory
        """
        
        # Create some folders
        fc_dirs = []
        for i in range(3):
            fc_dirs.append(tempfile.mkdtemp(suffix="{}_fc_dir".format(i), dir=self.rootdir))
        
        # Assert nothing is returned for non-existing flowcell
        self.assertListEqual(IlluminaRun._get_flowcell(self.rootdir,'flowcell_id'),
                             [],
                             "Nothing should be returned for a non-existing flowcell id")
        # Assert nothing is returned for None or empty flowcell id
        self.assertListEqual(IlluminaRun._get_flowcell(self.rootdir,''),
                             [],
                             "Nothing should be returned for a 0-length flowcell id")
        self.assertListEqual(IlluminaRun._get_flowcell(self.rootdir,None),
                             [],
                             "Nothing should be returned for an empty flowcell id")
        
        # Assert that the correct folders are returned for exact matches
        for fc_dir in fc_dirs:
            self.assertListEqual([fc_dir],
                                 IlluminaRun._get_flowcell(self.rootdir,os.path.basename(fc_dir)),
                                 "Did not return the correct folder for exact match")
        
        # Assert that a partial match is resolved to the correct folder
        self.assertListEqual([fc_dirs[-1]],
                             IlluminaRun._get_flowcell(self.rootdir,"{}_fc_dir".format(str(len(fc_dirs)-1))),
                             "Did not return the correct folder for partial match")
        
        # Assert that an ambiguous match returns the matching folders
        self.assertListEqual(sorted(fc_dirs),
                             sorted(IlluminaRun._get_flowcell(self.rootdir,"_fc_dir")),
                             "Did not return the correct folders for ambiguous matches")
        
        # Assert that the correct folder is returned for an exact match that matches ambiguously when allowing wildcards
        ambig_dir = os.path.join(self.rootdir,"_fc_dir")
        utils.safe_makedir(ambig_dir)
        self.assertListEqual([ambig_dir],
                             sorted(IlluminaRun._get_flowcell(self.rootdir,os.path.basename(ambig_dir))),
                             "Did not return the correct folder for specific non-wildcard match")
    
    def test__get_samplesheet(self):
        """Locate the samplesheet in a folder
        """
        # Create a few random files and folders and assert that they are not returned
        suffixes = [".csv","",""]
        for n in range(3):
            os.mkdir(os.path.join(self.rootdir,''.join(random.choice(string.ascii_uppercase) for x in range(5))))
            fh, _ = tempfile.mkstemp(dir=self.rootdir, suffix=suffixes[n])
            os.close(fh)
            
        self.assertIsNone(IlluminaRun._get_samplesheet(self.rootdir),
                          "Getting non-existing samplesheet did not return None")
        
        # Create a SampleSheet.csv and a [FCID].csv file and assert that they are
        # returned with a preference for the [FCID].csv file
        fcid = td.generate_fc_barcode()
        fcdir = os.path.join(self.rootdir,td.generate_run_id(fc_barcode=fcid))
        os.mkdir(fcdir)
        
        ss = [os.path.join(fcdir,"SampleSheet.csv"),
              os.path.join(fcdir,"{}.csv".format(fcid))]
        for s in ss:
            utils.touch_file(s)
            self.assertEqual(s,IlluminaRun._get_samplesheet(fcdir),
                             "Did not get existing {}".format(os.path.basename(s)))
    
