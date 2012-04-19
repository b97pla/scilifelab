
import os
import sys
import json
import operator
import time 
import subprocess
from optparse import OptionParser
import bcbio.distributed.slurm


def main(jobids, interval, log):
    
    # Continue monitoring until all the jobs have finished
    with open(log,"w") as logh:
    
        jobs_exist = True
        while jobs_exist:
            jobs_exist = False
            for jobid in jobids:
                jobs_exist |= poll_jobid(jobid, logh)
            time.sleep(interval)

def log_entry(entry, log):
    log.write("%s\n" % json.dumps(entry))
    
def poll_jobid(jobid, log):
    
    # Skip this jobid if it has finished
    if not bcbio.distributed.slurm.exists(jobid):
        return False
            
    # Get the status fields from squeue
    format = '%i;%j;%a;%P;%D;%C;%T;%N'
    status = _get_slurm_status(jobid,format)
    
    # If the job is running, get its resource consumption from the executing host
    cpu = [0,0]
    mem = [0,0]
    if status["STATE"] == "RUNNING":
        cpu, mem = _get_resource_usage(status["NODELIST"])
    (status["CPU_USER"],status["CPU_IDLE"]) = cpu  
    (status["MEM_USED"],status["MEM_FREE"]) = mem  
    
    entry = {time.time(): status}
    log_entry(entry,log)
    
    return True

def _get_resource_usage(host):
    
    # ssh to the host, execute top and parse the output
    cl = ["ssh",host]
    cl += ["top","-b","-n",2,"-d",5]
    cl += ["|","grep","-e","Cpu","-e","Mem"]
    cl = [str(i) for i in cl]
    output = subprocess.check_output(cl).split("\n")
    
    # Discard the first iteration and use the second for the resources
    cpu = {}
    mem = {}
    for line in output[2:]:
        if line.startswith("Cpu(s)"):
            parts = line[7:].split(",")
            for p in parts:
                c = p.strip().split('%')
                if len(c) != 2 or c[1] not in ("us","id"):
                    continue
                cpu[c[1]] = float(c[0])
        elif line.startswith("Mem"):
            parts = [p.strip() for p in line[4:].split(",")]
            for p in parts:
                m = p.split()
                if len(m) != 2 or m[1] not in ("used","free"):
                    continue
                mem[m[1]] = float(m[0][:-1])
    return [[cpu["us"],cpu["id"]],[mem["used"],mem["free"]]]
                
    
def _get_slurm_status(jobid, format):
    
    cl = ["squeue"]
    cl += ["-j",jobid]
    cl += ["-o","%s" % format]
    cl = [str(i) for i in cl]
    output = subprocess.check_output(cl).split("\n")
    
    # Create a dictionary from the first and second rows
    output = [line.split(";") for line in output[0:2]]
    status = dict(zip(output[0],output[1]))

    return status

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--interval", dest="interval", default=60)
    parser.add_option("-o", "--out_file", dest="out_file", default=sys.stdout)
    options, args = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()

    main(args,float(options.interval),options.out_file)
    