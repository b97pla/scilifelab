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
        self.batchsize = batchsize

    def __iter__(self):
        return self

    def next(self):
        l = min(len(self.filelist),self.batchsize)
        if l > 0:
            return [self.filelist.popleft() for i in range(l)]
        raise StopIteration

def run_screen(run_folder, batchsize, projectid, timelimit, jobname, email, slurm_extra):
    
    assert os.path.exists(run_folder), "The supplied run folder {} does not exist".format(run_folder)
    assert subprocess.check_call(["fastq_screen","-v"]) == 0, "fastq_screen could not be run properly, please check your environment"
    
    indirpattern = "*barcode"
    infilepattern = "*_1_fastq.txt*"
    
    # Setup outdir and bash script
    outdir = os.path.join(run_folder,"fastq_screen")
    try:
        os.mkdir(outdir,0770)
    except OSError:
        pass
    assert os.path.exists(outdir), "The output folder {} could not be created".format(outdir)
    
    bashscript = os.path.join(outdir,"run_fastq_screen.sh")
    with open(bashscript,"w") as fh:
        fh.write(run_script(outdir))
    os.chmod(bashscript,0770)
    
    sbatch_cmd = ["sbatch","--mail-user={}".format(email),
                  "--mail-type=FAIL","-D",outdir,"-A",projectid,"-J",jobname,"-t",timelimit,
                  "-N","1","-p","node"]
    sbatch_cmd += slurm_extra
    
    batch = Batcher(os.path.join(run_folder,indirpattern,infilepattern),batchsize)
    for infiles in batch:
        submit_batch(infiles,outdir,bashscript,sbatch_cmd)

def submit_batch(infiles, outdir, bashscript, sbatch_cmd):
    """Submits a slurm job that will run fastq_screen jobs on a common node 
       for the supplied infiles
    """ 
    sbatchfile = os.path.join(outdir,"fastq_screen_sbatch.sh")
    with open(sbatchfile,"w") as fh:
        fh.write("#! /bin/sh\n")
        for infile in infiles:
            pairfile = infile.replace("_1_fastq","_2_fastq")
            fh.write("{0} {1} {2} &\n".format(bashscript,infile,pairfile))
        fh.write("wait\n")
    os.chmod(sbatchfile,0770)
    
    sbatch_cmd += [sbatchfile]
    print subprocess.check_output(sbatch_cmd)
        
def run_script(outdir):
    """Return a bash script (as a text string) that will run fastq_screen"""
       
    return "#! /bin/sh\n" \
    "F1=$1\n" \
    "F2=$2\n" \
    "C=\"\"\n" \
    "if [ ${F1##*.} == \"gz\" ]\n" \
    "then\n" \
    "  C=\"compressed\"\n" \
    "  F1=${F1/.gz/}\n" \
    "  F2=${F2/.gz/}\n" \
    "  gzip -c -d ${F1}.gz > $F1 &\n" \
    "  if [ -e \"${F2}.gz\" ]\n" \
    "  then\n" \
    "    gzip -c -d ${F2}.gz > $F2 &\n" \
    "  fi\n" \
    "  wait\n" \
    "fi\n" \
    "fastq_screen --subset 2000000 --illumina --outdir %s --multilib $F1 --paired $F2\n" \
    "if [ \"$C\" == \"compressed\" ]\n" \
    "then\n" \
    "  if [ -e \"${F1}.gz\" ]\n" \
    "  then\n" \
    "    rm $F1\n" \
    "  fi\n" \
    "  if [ -e \"${F2}.gz\" ]\n" \
    "  then\n" \
    "    rm $F2\n" \
    "  fi\n" \
    "else\n" \
    "  gzip $F1 &\n" \
    "  gzip $F2 &\n" \
    "  wait\n" \
    "fi\n" % outdir
    
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
    parser.add_argument('-t','--time', action='store', default="0-06:00:00", 
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
    run_screen(args.run_folder[0], args.batchsize, args.projectid, args.time, args.jobname, args.email, args.slurm_extra.split())
      
if __name__ == "__main__":
    main()
