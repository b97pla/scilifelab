#!/usr/bin/env python
import drmaa
import os
import sys
import glob
import time
import yaml
import subprocess
from optparse import OptionParser

from bcbio.pipeline.alignment import get_genome_ref

REFERENCE_DIR = "/bubo/nobackup/uppnex/reference/biodata/galaxy"
SLURM_ARGS="-A a2010002 -p node"
RUN_TIME="6:00:00"

def main(input_path, genome, filter_file, read1, read2, filtered_reads, aligner, slurm_parameters):
    
    if filter_file is None:
        filter_file, _ = get_genome_ref(genome, aligner, os.path.normpath(REFERENCE_DIR))
    
    infiles = []
    if read1 is None:
        if os.path.isdir(input_path):
            pat = os.path.join(input_path,"*barcode","*_1_fastq.txt")
            for read1 in glob.glob(pat):
                read2 = read1.replace("_1_fastq.txt","_2_fastq.txt")
                if not os.path.exists(read2):
                    read2 = None
                infiles.append([read1,read2])
        elif os.path.isfile(input_path):
            if input_path.endswith("_1_fastq.txt"):
                read1 = input_path
                read2 = read1.replace("_1_fastq.txt","_2_fastq.txt")
                if not os.path.exists(read2):
                    read2 = None
            elif input_path.endswith("_2_fastq.txt"):
                read2 = input_path
                read1 = read1.replace("_2_fastq.txt","_1_fastq.txt")
                assert os.path.exists(read1), "ERROR: Could not find the first read file (expected %s)" % read1
            else:
                read1 = input_path
                read2 = None
            infiles.append([read1,read2])
    else:
        infiles.append([read1,read2])
    
    for read1, read2 in infiles:
        jobid = filter_files_job(read1, read2, filtered_reads,
                                 filter_file, aligner,
                                 slurm_parameters)
        print "Your job was submitted with jobid %s" % jobid

def filter_files_job(read1, read2, filtered_reads, filter_file, aligner, slurm_parameters):
    
    slurm_project, slurm_run_time, slurm_partition = slurm_parameters
    slurm_args = "-A %s -t %s -p %s" % (slurm_project,slurm_run_time,slurm_partition)
    
    print "Filtering %s against %s using %s" % (read1,filter_file,aligner)
    if read2 is not None:
        print "Including %s" % read2
    
    outname = read1.replace("fastq","-filtered.fastq")
    assert outname != read1, "ERROR: Could not properly interpret file name"
    
    args = ["--phred64-quals"]
    if slurm_partition == "node" or slurm_partition == "devel": args += ["-p", "8"]
    args += ["--un", outname,
             "-X", "6000"]
    if aligner == "bowtie2": args += ["-x"]
    args += [filter_file]
    if read2 is not None: 
        args += ["-1", read1,"-2", read2]
    else:
        args += [read1]
    if aligner == "bowtie2": args += ["-S"]
    if filtered_reads is not None:
        args += [filtered_reads]
    else:
        args += ["/dev/null"]
    
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = aligner   
    jt.args = args 
    
    jt.nativeSpecification = slurm_args

    jobid = s.runJob(jt)

    s.deleteJobTemplate(jt)
    s.exit()      
        
    return jobid
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-a", "--aligner", dest="aligner", default="bowtie")
    parser.add_option("-g", "--genome", dest="genome", default="phix")
    parser.add_option("-f", "--filter-file", dest="filter_file", default=None)
    parser.add_option("-1", "--read1", dest="read1", default=None)
    parser.add_option("-2", "--read2", dest="read2", default=None)
    parser.add_option("-r", "--filtered-reads", dest="filtered_reads", default=None)
    parser.add_option("-p", "--slurm-partition", dest="slurm_partition", default="core")
    parser.add_option("-A", "--slurm-project", dest="slurm_project", default="a2010002")
    parser.add_option("-t", "--slurm-run-time", dest="slurm_run_time", default="3:00:00")
    options, args = parser.parse_args()
    
    input_path = None
    if len(args) == 1:
        input_path = args[0]
    else:
        print __doc__
        sys.exit()
    main(os.path.normpath(input_path), 
         options.genome, options.filter_file, 
         options.read1, options.read2, options.filtered_reads,
         options.aligner, [options.slurm_project, 
                           options.slurm_run_time, 
                           options.slurm_partition])

