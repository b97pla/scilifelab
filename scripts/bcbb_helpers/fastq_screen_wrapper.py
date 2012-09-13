import glob
import os
import collections
import argparse
import subprocess

class Batcher:
    """Takes a glob and a batchsize and returns an iterator that will 
       return the resulting elements in batches of the specified size
    """
    def __init__(self, pattern, batchsize=4):
        self.filelist = collections.deque(glob.glob(pattern))
        assert batchsize > 0, "The batchsize needs to be > 0"
        self.batchsize = batchsize
        self._batchno = 0

    def __iter__(self):
        return self

    def next(self):
        l = min(len(self.filelist),self.batchsize)
        if l > 0:
            self._batchno += 1
            return [self.filelist.popleft() for i in range(l)]
        raise StopIteration
     
    def batchno(self):
        return self._batchno

def run_screen(run_folder, batchsize, projectid, timelimit, jobname, email, slurm_extra):
    
    assert os.path.exists(run_folder), "The supplied run folder {} does not exist".format(run_folder)
    assert subprocess.check_call(["fastq_screen","-v"]) == 0, "fastq_screen could not be run properly, please check your environment"
    
    indirpattern = "*/*XX/nophix"
    infilepattern = "*_1_fastq.txt*"
    
    outdir = run_folder
    bashscript = os.path.join(outdir,"run_fastq_screen.sh")
    with open(bashscript,"w") as fh:
        fh.write(run_script())
    os.chmod(bashscript,0770)
    
    sbatch_cmd = ["sbatch","--mail-user={}".format(email),
                  "--mail-type=FAIL","-D",outdir,"-A",projectid,"-J",jobname,"-t",timelimit,
                  "-N","1","-p","node"]
    sbatch_cmd += slurm_extra
    
    batch = Batcher(os.path.join(run_folder,indirpattern,infilepattern),batchsize)
    for infiles in batch:
        submit_batch(infiles,outdir,bashscript,sbatch_cmd,batch.batchno())

def submit_batch(infiles, outdir, bashscript, sbatch_cmd, batchno):
    """Submits a slurm job that will run fastq_screen jobs on a common node 
       for the supplied infiles
    """ 
    sbatchfile = os.path.join(outdir,"fastq_screen_sbatch_{}.sh".format(batchno))
    with open(sbatchfile,"w") as fh:
        fh.write("#! /bin/sh\n")
        for infile in infiles:
            pairfile = infile.replace("_1_fastq","_2_fastq")
            fh.write("{0} {1} {2} &\n".format(bashscript,infile,pairfile))
        fh.write("wait\n")
    os.chmod(sbatchfile,0770)
    
    print subprocess.check_output(sbatch_cmd + ["-o", "{}.out".format(os.path.splitext(os.path.basename(sbatchfile))[0]), sbatchfile]) 
        
def run_script():
    """Return a bash script (as a text string) that will run fastq_screen"""
       
    return "#! /bin/sh\n" \
    "F1=$1\n" \
    "F2=$2\n" \
    "OUTDIR=`dirname $F1`\"/../fastq_screen\"\n" \
    "mkdir -p $OUTDIR\n" \
    "if [ ${F1##*.} == \"gz\" ]\n" \
    "then\n" \
    "  F1=${F1/.gz/}\n" \
    "  gzip -c -d ${F1}.gz > $F1 &\n" \
    "fi\n" \
    "if [ ${F2##*.} == \"gz\" ]\n" \
    "then\n" \
    "  F2=${F2/.gz/}\n" \
    "  gzip -c -d ${F2}.gz > $F2 &\n" \
    "fi\n" \
    "wait\n" \
    "fastq_screen --subset 2000000 --outdir $OUTDIR --multilib $F1 --paired $F2\n" \
    "if [ -e ${F1}.gz ]\n" \
    "then\n" \
    "  rm $F1\n" \
    "else\n" \
    "  gzip $F1 &\n" \
    "fi\n" \
    "if [ -e ${F2}.gz ]\n" \
    "then\n" \
    "  rm $F2\n" \
    "else\n" \
    "  gzip $F2 &\n" \
    "fi\n" \
    "wait\n"
    
def main():
    
    parser = argparse.ArgumentParser(description="Run fastq_screen on the demultiplexed fastq "\
                                     "files in the *_barcode directories in a run folder. "\
                                     "The jobs will be batched together in groups of [batchsize] "\
                                     "(default=4) and submitted to slurm as a node job. If the "\
                                     "input files are gzipped, a temporary decompressed copy "\
                                     "is placed in the *_barcode directory and removed after "\
                                     "processing. If the input files are uncompressed, they will "\
                                     "be gzipped after processing.")

    parser.add_argument('-b','--batchsize', action='store', default=4, 
                        help="the number of fastq_screen jobs to run together on one node")
    parser.add_argument('-A','--projectid', action='store', default="a2010002", 
                        help="the Uppnex project id to use with slurm")
    parser.add_argument('--time', action='store', default="0-06:00:00", 
                        help="the time limit for the slurm job")
    parser.add_argument('-J','--jobname', action='store', default="fastq_screen", 
                        help="the job name to use for slurm jobs")
    parser.add_argument('-e','--email', action='store', default="seqmaster@scilifelab.se", 
                        help="email address for slurm reporting")
    parser.add_argument('--slurm_extra', action='store', default="", 
                        help="extra parameters to pass to sbatch")
    parser.add_argument('run_folder', action='store', default=None, 
                        help="the full path to the run folder containing analysis output", nargs='+')
    
    args = parser.parse_args()
    run_screen(os.path.abspath(args.run_folder[0]), args.batchsize, args.projectid, args.time, args.jobname, args.email, args.slurm_extra.split())
      
if __name__ == "__main__":
    main()


# --- Tests --- #

import unittest
import tempfile
import shutil
from mock import Mock

class TestBatcher(unittest.TestCase):
    
    def setUp(self):
        """Set up a test folder for Batcher
        """
        self._test_folder = tempfile.mkdtemp(prefix="BatcherTest")
        
    def tearDown(self):
        shutil.rmtree(self._test_folder)
        
    def test_1(self):
        """TestBatcher: Test a setsize smaller than batchsize
        """
        setsize, batchsize = 3,4
        self._test_batch(setsize,batchsize)
        
    def test_2(self):
        """TestBatcher: Test a setsize larger than batchsize
        """
        setsize, batchsize = 3,2
        self._test_batch(setsize,batchsize)
        
    def test_3(self):
        """TestBatcher: Test a setsize equal to batchsize
        """
        setsize, batchsize = 5,5
        self._test_batch(setsize,batchsize)
        
    def test_4(self):
        """TestBatcher: Test a setsize equal to 0
        """
        setsize, batchsize = 0,10
        self._test_batch(setsize,batchsize)
        
    def test_5(self):
        """TestBatcher: Test a setsize equal to twice the batchsize
        """
        setsize, batchsize = 10,5
        self._test_batch(setsize,batchsize)
        
    def test_6(self):
        """TestBatcher: Test a batchsize equal to 0
        """
        setsize, batchsize = 10,0
        with self.assertRaises(AssertionError):
            self._test_batch(setsize,batchsize)
            
    def _test_batch(self, setsize, batchsize):
        # Create a set of testfiles
        dummyfiles = []
        for i in range(setsize):
            fh, dummyfile = tempfile.mkstemp(suffix=".batcher_test", dir=self._test_folder)
            os.close(fh)
            dummyfiles.append(dummyfile)
        # Create a file that should not be returned
        fh, nomatchfile = tempfile.mkstemp(suffix=".batcher_nomatch_test", dir=self._test_folder)
        os.close(fh)
        
        # Pattern 
        pattern = os.path.join(self._test_folder, "*.batcher_test")
        batcher = Batcher(pattern,batchsize)
        
        expected_bno = 0
        for b in batcher:
            expected_size = min(batchsize, setsize - batchsize*expected_bno)
            expected_bno += 1
            # Assert that each batch has the expected size
            self.assertEqual(len(b),expected_size, "The returned batchsize does not match the expected size ({} != {})".format(len(b),expected_size))
            # Assert that the batch enumerator gives the correct batch number
            self.assertEqual(batcher.batchno(),expected_bno, "The returned batch number does not match the expected ({} != {})".format(batcher.batchno(),expected_bno))
    
            for bf in b:
                # Assert that the returned files existed in the original array and has not already been returned 
                self.assertIn(bf, dummyfiles, "The returned file is not present in the expected set ({})".format(bf))
                self.assertTrue(os.path.exists(bf), "The returned file has already been seen ({})".format(bf))
                # Remove the returned files
                os.unlink(bf)
        
        # Assert that all files have been returned by the iterator
        for dummyfile in dummyfiles:
            self.assertFalse(os.path.exists(dummyfile), "Some files were not iterated over ({})".format(dummyfile))
        os.unlink(nomatchfile)
        
class TestRunScreen(unittest.TestCase):
    
    def setUp(self):
        self.testdir = tempfile.mkdtemp(prefix="run_screen_test")
    
    def tearDown(self):
        shutil.rmtree(self.testdir)
        
    def test_missing_run_folder(self):
        """TestRunScreen: Test a missing run folder
        """
        run_folder = os.path.join(self.testdir,"non-existing-folder")
        with self.assertRaises(AssertionError):
            run_screen(run_folder, 0, "", "", "", "", "")
            
    def test_failing_fastq_screen(self):
        """TestRunScreen: Test that failing call to fastq_screen raises error
        """
        subprocess.check_call = Mock(return_value=1)
        with self.assertRaises(AssertionError):
            run_screen(self.testdir, 5, "", "", "", "", "")
            
    def test_bash_script_creation(self):
        """TestRunScreen: Test that bash script for running fastq_screen is created
        """
        subprocess.check_call = Mock(return_value=0)
        run_screen(self.testdir, 5, "", "", "", "", "")
        bashscript = os.path.join(self.testdir,"run_fastq_screen.sh")
        self.assertTrue(os.path.exists(bashscript))
        
    def test_sbatch_submission(self):
        """TestRunScreen: Test that call to sbatch submission is correct
        """
        subprocess.check_call = Mock(return_value=0)
        
        folders = os.path.join(self.testdir,"1","2CXX","nophix")
        os.makedirs(folders)
        for i in range(2):
            f = os.path.join(folders,"{}_1_fastq.txt.gz".format(str(i)))
            open(f,"w").close()
        
        # FIXME: Create mock object for submit_batch and check calling parameters
        #run_screen(self.testdir, 5, "", "", "", "", "")
        
        
        