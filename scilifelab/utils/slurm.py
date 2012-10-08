"""Useful functions for interacting with the slurm manager
"""

import subprocess
import getpass
import drmaa

def get_slurm_jobid(jobname,user=getpass.getuser()):
    """Attempt to get the job id for a slurm job name. Can this be done with python-drmaa instead?
    """
    jobids = []
    cmd = ['/usr/bin/squeue','-h','-o','%i','-n',jobname,'-u',user]
    try:
        retval = str(subprocess.check_output(cmd))
        for val in retval.split("\n"):
            jobids.append(int(val))
    except:
        pass
    return jobids

def get_slurm_jobstatus(jobid):
    """Get the status for a jobid
    """
    s = drmaa.Session()
    s.initialize()
    status = s.jobStatus(str(jobid))
    s.exit()
    return status
    
