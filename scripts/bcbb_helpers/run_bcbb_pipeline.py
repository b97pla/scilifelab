#!/usr/bin/env python
import drmaa
import os
import sys
import glob
import time
import yaml
import subprocess

ANALYSIS_SCRIPT="distributed_nextgen_pipeline.py"
ANALYSIS_DIR="/proj/a2010002/nobackup/illumina/"
EMAIL="seqmaster@scilifelab.se"
SLURM_ARGS="-A a2010002 -p node --qos=seqver"
RUN_TIME="168:00:00"
PROCESS_YAML_SCRIPT = "process_run_info.py"
PROCESS_YAML = True
FC_SPECIFIC_AMPQ = True

def main():
 
    post_process_config = sys.argv[1]
    archive_dir = sys.argv[2]
    if len(sys.argv) > 3:
        run_info_file = sys.argv[3]
    else:
        run_info_file = None
    
    # Set the barcode type in run_info.yaml to "illumina", strip the 7th nucleotide and set analysis to 'Minimal'
    if run_info_file is not None and PROCESS_YAML:
        print "---------\nProcessing run_info:"
        run_info_backup = "%s.orig" % run_info_file
        os.rename(run_info_file,run_info_backup)
        cl = ["python","%s" % PROCESS_YAML_SCRIPT,run_info_backup,"--analysis","Align_illumina","--out_file",run_info_file,"--ascii","--clear_description"]
        print subprocess.check_output(cl)
        print "\n---------\n"
    
    # Check that the specified paths exist
    print "Checking input paths"
    for path in (post_process_config,archive_dir,run_info_file):
        if path is not None and not os.path.exists(path):
            raise Exception("The path %s does not exist" % path)
 
    print "Getting base_dir from %s" % post_process_config
    # Parse the config to get the analysis directory
    with open(post_process_config) as ppc:
        config = yaml.load(ppc)
    
    analysis = config.get("analysis",{})
    base_dir = analysis.get("base_dir",ANALYSIS_DIR)
    
    print "Getting run name from %s" % archive_dir
    # Get the run name from the archive dir
    _,run_name = os.path.split(os.path.normpath(archive_dir))

    # Create the working directory if necessary and change into it
    work_dir = os.path.join(base_dir,run_name)
    os.chdir(base_dir)
    print "Creating/changing to %s" % work_dir
    try:
        os.mkdir(run_name,0770)
    except OSError:
        pass
    os.chdir(run_name)
 
    # make sure that the work dir exists
    if not os.path.exists(work_dir):
        raise Exception("The path %s does not exist and was not created" % work_dir)
    
    # if required, parse the machine id and flowcell position and use an ampq vhost specific for it
    if FC_SPECIFIC_AMPQ:
        machine_id = None
        flowcell_position = None
        for p in run_name.upper().split("_"):
            if p.startswith("SN"):
                machine_id = p
            elif p[0] in ("A","B") and p.endswith("XX"):
                flowcell_position = p[0]
        assert machine_id and flowcell_position, "Machine id and flowcell position could not be parsed from run name '%s'" % run_name
        
        # write a dedicated post_process.yaml for the ampq queue
        if config.get("distributed",False):
            config["distributed"]["rabbitmq_vhost"] = "bionextgen-%s-%s" % (machine_id,flowcell_position)
        
        post_process_config_orig = post_process_config
        parts = os.path.splitext(post_process_config)
        post_process_config = "%s-%s-%s%s" % (parts[0],machine_id,flowcell_position,parts[1])
        
        if os.path.exists(post_process_config):
            print "Customized post process yaml '%s' already exists, will not overwrite" % post_process_config
        else:
            with open(post_process_config,"w") as fh:
                fh.write(yaml.safe_dump(config, default_flow_style=False, allow_unicode=True)) 
        
    # Translate the time string into seconds
    multiplier = (86400,3600,60,1)
    units = RUN_TIME.split(":")
    seconds = 0
    for i in range(1,len(units)+1):
        seconds += multiplier[-1*i]*int(units[-1*i])
    print "Will request job to run for %s seconds" % seconds
    print "Initializing session"
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = ANALYSIS_SCRIPT
    args = [post_process_config,archive_dir]
    if run_info_file is not None:
        args.append(run_info_file)
    jt.args = args
    
    # TODO: job name is always (null), must fix slurm_drmaa C library and its
    # custom parsing (substitute "slurmdrmaa_parse_native"
    # for GNU GetOpt on slurm_drmaa/util.c)
    #jt.job_name = run_name
    #jt.blockEmail = 0
    #jt.mail = EMAIL
    #jt.hardRunDurationLimit(seconds)
    output = os.path.join(work_dir,"%s.out" % run_name)
    error = os.path.join(work_dir,"%s.err" % run_name)
    print "Setting output and error paths (%s,%s)" % (output,error)
    
    #jt.outputPath = output
    #jt.errorPath = error 
    jt.nativeSpecification = "%s -t %s" % (SLURM_ARGS,RUN_TIME)

    print "Submitting job"
    jobid = s.runJob(jt)
    print 'Your job has been submitted with id ' + jobid

    s.deleteJobTemplate(jt)
    s.exit()

    exit(0)

if __name__ == "__main__":
    main()
