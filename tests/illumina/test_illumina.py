
import unittest
import shutil
import tempfile
import os
import random
import string
import bcbio.utils as utils
import tests.generate_test_data as td
import scilifelab.illumina as illumina
from scilifelab.illumina import IlluminaRun
from scilifelab.illumina import map_index_name
        
class TestIlluminaRun(unittest.TestCase):
    
    def setUp(self): 
        self.rootdir = tempfile.mkdtemp(prefix="test_illumina_run_")
        
        # Create a fcdir
        self.exp_fcid = td.generate_fc_barcode()
        self.exp_fcdir = os.path.join(self.rootdir,td.generate_run_id(fc_barcode=self.exp_fcid))
        os.mkdir(self.exp_fcdir)
        
        # Create multiple sequence read directories
        self.exp_seqdir = [os.path.join(self.exp_fcdir,"Unaligned"),
                           os.path.join(self.exp_fcdir,"Unaligned_L6"),
                           os.path.join(self.exp_fcdir,"Unaligned_L8")]
        for d in self.exp_seqdir:
            os.mkdir(d)
        
        # Create directories for undetermined indices reads
        self.exp_unmatched_directory = [os.path.join(d,"Undetermined_indices") for d in self.exp_seqdir[:-1]]
        for d in self.exp_unmatched_directory:
            os.mkdir(d)

        # Create basecall stats directories
        self.exp_basecall_stats = [os.path.join(d,"Basecall_Stats_{}".format(self.exp_fcid)) for d in self.exp_seqdir[1:]]
        for d in self.exp_basecall_stats:
            os.mkdir(d)
            
        self.run = IlluminaRun(self.exp_fcdir)
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)

    def test_map_index_name(self):
        """Map index sequences to names
        """
        
        # Some test cases
        from scilifelab.illumina.index_definitions import BASIC_LOOKUP
        random_keys = [BASIC_LOOKUP.keys()[random.randint(0,len(BASIC_LOOKUP.keys())-1)] for n in xrange(25)]
        
        for name in random_keys:
            self.assertIn(name,map_index_name(BASIC_LOOKUP[name],0),
                          "Exact mapping did not return expected index name")
      
    def test_get_flowcell(self):
        """Get flowcell from analysis directory
        """
        
        # Work in a separate subdirectory
        subdir = os.path.join(self.rootdir,"test_get_flowcell")
        os.mkdir(subdir)
        
        # Create some folders
        fc_dirs = []
        for i in range(3):
            fc_dirs.append(tempfile.mkdtemp(suffix="{}_fc_dir".format(i), dir=subdir))
        
        # Assert nothing is returned for non-existing flowcell
        self.assertListEqual(IlluminaRun.get_flowcell(subdir,'flowcell_id'),
                             [],
                             "Nothing should be returned for a non-existing flowcell id")
        
        # Assert that the correct folders are returned for exact matches
        for fc_dir in fc_dirs:
            self.assertListEqual([fc_dir],
                                 IlluminaRun.get_flowcell(subdir,os.path.basename(fc_dir)),
                                 "Did not return the correct folder for exact match")
        
        # Assert that an empty match returns all folders
        self.assertListEqual(sorted(fc_dirs),
                             sorted(IlluminaRun.get_flowcell(subdir)),
                             "Did not return the correct folders for empty matches")
        
        # Assert that a partial match is resolved to the correct folder
        self.assertListEqual([fc_dirs[-1]],
                             IlluminaRun.get_flowcell(subdir,"{}_fc_dir".format(str(len(fc_dirs)-1))),
                             "Did not return the correct folder for partial match")
        
        # Assert that an ambiguous match returns the matching folders
        self.assertListEqual(sorted(fc_dirs),
                             sorted(IlluminaRun.get_flowcell(subdir,"_fc_dir")),
                             "Did not return the correct folders for ambiguous matches")
        
        # Assert that the correct folder is returned for an exact match that matches ambiguously when allowing wildcards
        ambig_dir = os.path.join(subdir,"_fc_dir")
        utils.safe_makedir(ambig_dir)
        self.assertListEqual([ambig_dir],
                             sorted(IlluminaRun.get_flowcell(subdir,os.path.basename(ambig_dir))),
                             "Did not return the correct folder for specific non-wildcard match")
        
        shutil.rmtree(subdir)
    
    def test_get_samplesheet(self):
        """Locate the samplesheet in a folder
        """
        
        # Work in a separate subdirectory
        subdir = os.path.join(self.rootdir,"test_get_samplesheet")
        os.mkdir(subdir)
        
        # Create a few random files and folders and assert that they are not returned
        suffixes = [".csv","",""]
        for n in range(3):
            os.mkdir(os.path.join(subdir,''.join(random.choice(string.ascii_uppercase) for x in range(5))))
            fh, _ = tempfile.mkstemp(dir=subdir, suffix=suffixes[n])
            os.close(fh)
            
        self.assertIsNone(IlluminaRun.get_samplesheet(subdir),
                          "Getting non-existing samplesheet did not return None")
        
        # Create a SampleSheet.csv and a [FCID].csv file and assert that they are
        # returned with a preference for the [FCID].csv file
        fcid = td.generate_fc_barcode()
        fcdir = os.path.join(subdir,td.generate_run_id(fc_barcode=fcid))
        os.mkdir(fcdir)
        
        ss = [os.path.join(fcdir,"SampleSheet.csv"),
              os.path.join(fcdir,"{}.csv".format(fcid))]
        for s in ss:
            utils.touch_file(s)
            self.assertEqual(s,IlluminaRun.get_samplesheet(fcdir),
                             "Did not get existing {}".format(os.path.basename(s)))
    
        shutil.rmtree(subdir)
        
    def test_get_sequence_dir(self):
        """Get the sequence read directory
        """
        self.assertEqual(self.exp_seqdir, self.run.get_sequence_dir(),
                         "Did not get correct top directory for sequence reads")
        
    def test_get_unmatched_dir(self):
        """Get the unmatched indices read top directory
        """
        self.assertEqual(self.exp_unmatched_directory, self.run.get_unmatched_dir(),
                         "Did not get correct top directory for undetermined indices sequence reads")
        
    def test_get_basecall_stats(self):
        """Get the basecall stats directory
        """
        self.assertEqual(self.run.get_basecall_stats(), self.exp_basecall_stats,
                         "Did not get correct Basecall_Stats directory")
        
    def test_get_unmatched_reads(self):
        """Get the undetermined indexes reads
        """
        
        # Create some files representing undetermined index reads
        lanes = [1,3,5,7]
        readfiles = []
        for lane in lanes:
            fdir = os.path.join(self.exp_unmatched_directory[0],"Sample_lane{:d}".format(lane))
            readfiles.append([os.path.join(fdir,"lane{l:d}_Undetermined_L00{l:d}_R{r:d}_*.fastq.gz".format(l=lane, r=read)) for read in [1,2]])
            os.makedirs(fdir)
            for readfile in readfiles[-1]:
                utils.touch_file(readfile)
        
            # Assert that the correct files are returned
            self.assertListEqual(sorted(readfiles[-1]),
                                 sorted(self.run.get_unmatched_reads(lanes=[lane])[0]),
                                 "Did not get expected undetermined indexes reads")

             