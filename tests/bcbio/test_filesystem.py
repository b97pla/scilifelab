
import os
import subprocess
import unittest
import tempfile
import shutil
import random
import string
import datetime
import bcbio.utils as utils
import tests.generate_test_data as td
import scilifelab.bcbio.filesystem as sq

class TestFilesystem(unittest.TestCase):
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_bcbio_filesystem_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
        
    def test__get_project_analysis_dir(self):
        """Test that getting the project analysis folder behaves as expected
        """
        # Assert that none is returned when no folder exists
        proj = td.generate_project()
        self.assertIsNone(sq.get_project_analysis_dir(self.rootdir,proj),
                          "Did not return empty result for non-existing folders")
        
        # Assert that none is still returned when some mismatching folders exist
        for n in range(5):
            os.mkdir(os.path.join(self.rootdir,td.generate_project()))
        self.assertIsNone(sq.get_project_analysis_dir(self.rootdir,proj),
                          "Did not return empty result for mismatching folders")
        
        # Assert that a file with the same name as the project is not returned
        projdir = os.path.join(self.rootdir,proj)
        utils.touch_file(projdir)
        self.assertIsNone(sq.get_project_analysis_dir(self.rootdir,proj),
                          "Returned a file with matching name. Should only return folders")
        os.unlink(projdir)
        
        # Assert that the corrct folder is returned when it exists
        os.mkdir(projdir)
        self.assertEqual(projdir,sq.get_project_analysis_dir(self.rootdir,proj),
                         "The expected project folder was not returned")
        
    def test_get_sample_pipeline_log(self):
        """Get the sample pipeline logfile
        """
        sample = td.generate_sample()
        # Assert that an empty directory returns None
        self.assertIsNone(sq.get_sample_pipeline_log(self.rootdir,sample),
                          "Getting a sample log from an empty directory did not return None")
        
        # Assert that non-relevant log files are not returned
        for n in range(5):
            utils.touch_file(os.path.join(self.rootdir,"{}-bcbb.log".format(td.generate_sample())))
        self.assertIsNone(sq.get_sample_pipeline_log(self.rootdir,sample),
                          "Getting a non-existing sample log did not return None")
        
        # Assert that the correct log file is returned when it exists
        utils.touch_file(os.path.join(self.rootdir,"{}-bcbb.log".format(sample)))
        self.assertEqual(os.path.join(self.rootdir,"{}-bcbb.log".format(sample)),sq.get_sample_pipeline_log(self.rootdir,sample),
                         "Getting an existing sample log file did not return the expected file")
        
    def test_get_pipeline_indicator(self):
        """Get pipeline indicator files
        """
        
        self.assertListEqual([],sq.get_pipeline_indicator(self.rootdir),
                          "Empty directory did not return an empty list")
        
        # Create some random files and assert that they are not picked up
        for n in range(5):
            fh, _ = tempfile.mkstemp(dir=self.rootdir)
            os.close(fh)
        self.assertListEqual([],sq.get_pipeline_indicator(self.rootdir),
                          "Non-existing indicator files did not return an empty list")
        
        # Create some indicator files and assert that they are returned
        ifiles = []
        for n in range(1,6):
            ifiles.append(os.path.join(self.rootdir,"{s:02d}_{act}.txt".format(s=n,act=''.join(random.choice(string.ascii_lowercase) for x in range(5)))))
            utils.touch_file(ifiles[-1])
            
        self.assertListEqual(sorted(ifiles),sorted(sq.get_pipeline_indicator(self.rootdir)),
                          "Existing indicator files did not return the expected output")
        
        # Assert that asking for a specific indicator returns the expected output
        self.assertListEqual(sorted(ifiles[0:2]),sorted(sq.get_pipeline_indicator(self.rootdir,range(1,3))),
                          "Specific indicator files did not return the expected output")
        
        # Assert that asking for a specific non-existing indicator returns the expected output
        self.assertListEqual([ifiles[-1]],sorted(sq.get_pipeline_indicator(self.rootdir,range(len(ifiles),len(ifiles)+2))),
                          "Specific non-existing indicator files did not return the expected output")
        
    def test_get_most_recent_indicator(self):
        """Get the most recent timestamp and corresponding file
        """
        
        # Assert that an empty input returns a (0,None) tuple
        self.assertEqual((datetime.datetime.fromtimestamp(0.0),None),sq.get_most_recent_indicator([]),
                             "No input files did not return expected output")
        
        # Assert that a file without timestamp returns (0,None)
        fh, tfile = tempfile.mkstemp(dir=self.rootdir)
        os.close(fh)
        self.assertEqual((datetime.datetime.fromtimestamp(0.0),None),sq.get_most_recent_indicator([tfile]),
                             "Empty input file did not return expected output")
        
        # Assert that the file with the most recent timestamp is returned
        tfiles = []
        for t in [1000.,2000.,3000.,4000.]:
             fh, tfile = tempfile.mkstemp(dir=self.rootdir)
             os.write(fh,"{}\n".format(datetime.datetime.fromtimestamp(t).isoformat()))
             os.close(fh)
             tfiles.append(tfile)
        
        self.assertEqual((datetime.datetime.fromtimestamp(4000.0),tfiles[-1]),sq.get_most_recent_indicator(tfiles),
                             "Timestamped input files did not return expected output")
        
        with open(tfiles[0],"a") as fh:
            fh.write("{}\n".format(datetime.datetime.fromtimestamp(5000.0).isoformat()))
        self.assertEqual((datetime.datetime.fromtimestamp(5000.0),tfiles[0]),sq.get_most_recent_indicator(tfiles),
                             "Input files with multiple timestamps did not return expected output")
            
    def test_fastq_screen_finished(self):
        """Detecting finished state of fastq_screen
        """
        
        # Assert that an empty directory doesn't indicate finished state
        self.assertFalse(sq.fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished without output files")
        
        # Create an output file and corresponding png but no rows in output
        sample_file = os.path.join(self.rootdir,"{}_fastq_screen.txt".format(td.generate_sample()))
        png_file = "{}.png".format(os.path.splitext(sample_file)[0])
        utils.touch_file(sample_file)
        utils.touch_file(png_file)
        self.assertFalse(sq.fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished with empty output file")
        
        # Write some output and assert fastq_screen is detected as finished
        with open(sample_file,"w") as fh:
            for n in range(5):
                fh.write("{}\n".format(str(n)))
        self.assertTrue(sq.fastq_screen_finished(self.rootdir),
                         "Fastq screen should be considered finished with non-empty output file and corresponding png")
        
        # Remove the png and assert fastq_screen is not finished
        os.unlink(png_file)
        self.assertFalse(sq.fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished with non-empty output file but without corresponding png")
        
