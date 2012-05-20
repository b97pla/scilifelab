#!/usr/bin/env python
"""
    run_cmd.py cmd args [slurm options]
"""
import drmaa
import os
import sys
import glob
import time
import subprocess
from optparse import OptionParser

def main(cmd, args, slurm_parameters):
    slurm_project, slurm_run_time, slurm_partition, slurm_job_name, slurm_output, slurm_error, slurm_opts = slurm_parameters
    slurm_args = "-A %s -t %s -p %s -J %s" % (slurm_project,slurm_run_time,slurm_partition,slurm_job_name )
    # if not slurm_output is None:
    #     slurm_args += " -o %s" % slurm_output
    # if not slurm_error is None:
    #     slurm_args += " -e %s" % slurm_error
    if not slurm_opts is None:
        slurm_args += " %s" % slurm_opts
        
    print "Command: %s" % cmd
    print "Args: %s" % args
    print "slurm args: %s" % slurm_args
    
    if options.dry_run:
        sys.exit()
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = cmd
    jt.args = args
    
    jt.nativeSpecification = slurm_args

    jobid = s.runJob(jt)

    s.deleteJobTemplate(jt)
    s.exit()

    print "Your job was submitted with jobid %s" % jobid    
    
    
if __name__ == "__main__":
    usage = """
    run_cmd.py cmd args [slurm options]
    """
    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--slurm-partition", dest="slurm_partition", default="core")
    parser.add_option("-A", "--slurm-project", dest="slurm_project", default="a2010002")
    parser.add_option("-t", "--slurm-run-time", dest="slurm_run_time", default="3:00:00")
    parser.add_option("-J", "--slurm-job-name", dest="slurm_job_name", default="run_cmd")
    parser.add_option("-e", "--slurm-error", dest="slurm_error", default=None)
    parser.add_option("-o", "--slurm-output", dest="slurm_output", default=None)
    parser.add_option("-O", "--slurm-opts", dest="slurm_opts", default=None)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, action="store_true")
    parser.add_option("-n", "--dry_run", dest="dry_run", default=False, action="store_true")
    options, args = parser.parse_args()

    print "Number of arguments: %s" % len(args)    
    if len(args) != 1:
        print __doc__
        sys.exit()

    kwargs = dict(
        cmd = args[0].split()[0],
        args=args[0].split()[1:len(args[0].split())],
        slurm_parameters =[options.slurm_project,
                           options.slurm_run_time,
                           options.slurm_partition,
                           options.slurm_job_name,
                           options.slurm_output,
                           options.slurm_error,
                           options.slurm_opts,
                           ],
        )
        
    main(**kwargs) 
