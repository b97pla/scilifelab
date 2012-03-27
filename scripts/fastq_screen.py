#!/usr/bin/env python
import drmaa
import os
import sys
import glob
import time

def main():

    samples = glob.glob(sys.argv[1]+"/*_barcode/*1_fastq.txt")

    for sample1 in samples:
        sample2 = sample1.replace("1_fastq", "2_fastq")

        s = drmaa.Session()
        s.initialize()

        jt = s.createJobTemplate()
        jt.remoteCommand = 'fastq_screen'
        jt.args = ['--illumina', '--multilib', '--subset', '2000000', sample1, \
                   '--paired', sample2, '--outdir', '.']

        # TODO: job name is always (null), must fix slurm_drmaa C library and its
        # custom parsing (substitute "slurmdrmaa_parse_native"
        # for GNU GetOpt on slurm_drmaa/util.c)
        jt.job_name = "fastq_screen %s" % os.path.basename(sample1)
        jt.nativeSpecification = "-A a2010002 -p node -t 00:30:00"

        jobid = s.runJob(jt)
        print 'Your job has been submitted with id ' + jobid

        s.deleteJobTemplate(jt)
        s.exit()

    exit(0)

if __name__ == "__main__":
    main()
