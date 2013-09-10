######
# pulse_conversion.py
#
# Written by: Pontus Larsson, Science for Life Laboratory, Stockholm
# 
# This is a proof-of-concept script performing bcl-to-fastq conversion 
# during the progress of an illumina run using the CASAVA software suite.
# It is not to be considered stable so use it at your own risk.
# 
# Usage: 
#       python pulse_conversion.py ILLUMINA_FLOWCELL_DIR [remote_user_name] [remote_host] [remote_path]
#
# If remote credentials are supplied, the resulting fastq files will be uploaded
# (assuming the SSH-keys have been set up correctly). Note that the CASAVA wrapper
# bash script 'convertBclToFastq.sh' (available in this GitHub repository) is required.
#
######

import xml.etree.ElementTree as ET
import sys
import os
import subprocess
import datetime
import smtplib
from email.mime.text import MIMEText

#### The parameters in this section can (and should) be customized
extraction_points = [35,50,75,100]
bclToFastq = "convertBclToFastq.sh"
cores = "8"
recipients = ["first_recipient@scilifelab.se","second_recipient@scilifelab.se"]
sender = "sender@scilifelab.se"
####

def _get_cycle(xmlfile):
    """Parse the config.xml and return the last cycle 
    """
    cycle = None
    if os.path.exists(xmlfile):
        tree = ET.ElementTree()
        tree.parse(xmlfile)
        element = tree.find("Run/Cycles")
        cycle = int(element.get("Last"))
    return cycle

def _get_run_info(xmlfile):
    """Parse the RunInfo.xml and return the read setup, flowcell id and lane count
    """
    reads = []
    fcid = None
    lanes = 0
    if os.path.exists(xmlfile):
        tree = ET.ElementTree()
        tree.parse(xmlfile)
        element = tree.find("Run/Reads")
        for read in element:
            reads.append(dict(zip(read.keys(), [read.get(k) for k in read.keys()])))
        element = tree.find("Run/Flowcell")
        fcid = element.text
        element = tree.find("Run/FlowcellLayout")
        lanes = int(element.get("LaneCount","0"))
        
    return sorted(reads, key=lambda r: int(r.get("Number", 0))), fcid, lanes

def _generate_samplesheet(samplesheet, fcid, lanes, cycle):
    header = ["FCID","Lane","SampleID","SampleRef","Index","Description","Control","Recipe","Operator","SampleProject"]
    with open(samplesheet,"w") as fh:
        fh.write("{}\n".format(",".join(header)))
        for lane in xrange(lanes):
            fh.write("{}\n".format(",".join([fcid,str(lane+1),"Cycle{}".format(str(cycle)),"","","","","","","pulse_analysis"])))

starttime = datetime.datetime.utcnow()
output = []

rundir = sys.argv[1]
destination = None
if len(sys.argv) >= 5:
    destination = "{}@{}:{}".format(*sys.argv[2:5])

# Change to the rundir
os.chdir(rundir)

# check that the config file exists
xmlfile = os.path.join(rundir,"RunInfo.xml")
assert os.path.exists(xmlfile), "Required file {} does not exist".format(xmlfile)

# get run setup
reads, fcid, lanes = _get_run_info(xmlfile)
assert len(reads) > 0, "Could not determine run setup"
assert fcid is not None, "Could not determine flowcell id"
assert lanes > 0, "Could not determine number of lanes on flowcell"

# second read offset
offset = sum([int(reads[i].get("NumCycles","0")) for i in xrange(len(reads)-1)])

# check that the config file exists
xmlfile = os.path.join(rundir,"Data","Intensities","BaseCalls","config.xml")
if not os.path.exists(xmlfile):
    #, "Required file {} does not exist".format(xmlfile)
    sys.exit(0)

# get the last cycle basecalled
cycle = _get_cycle(xmlfile)
assert cycle is not None, "Could not determine last cycle basecalled"

use_cycle = None
for ep in sorted([e + offset for e in extraction_points] + extraction_points, reverse=True):

    # Skip if we haven't reached the desired cycle
    if cycle <= ep:
        continue 

    # Break at the first (=the highest) extraction point that can be processed
    use_cycle = ep
    break
    
epfile = os.path.join(rundir,"Cycle{}_extracted.txt".format(use_cycle))
# Only proceed if this cycle hasn't already been extracted.
if use_cycle is not None and not os.path.exists(epfile):
    
    # Write the indicator file to indicate processing started
    open(epfile,"w").close()
    
    output.append("{}\tProcessing run {} at cycle {}".format(datetime.datetime.utcnow().isoformat(),rundir,use_cycle))
      
    # generate samplesheet
    samplesheet = os.path.join(rundir,"{}_pulse.csv".format(fcid))
    _generate_samplesheet(samplesheet,fcid,lanes,use_cycle)
    
    # basemask to use for second read
    basemask = ",".join(["N*" for i in xrange(len(reads)-1)])
    if use_cycle > int(reads[0].get("NumCycles","0")):
        basemask = "{},{}".format(basemask,"Y{}N*".format(str(use_cycle-offset)))
    else:
        basemask = "{},{}".format("Y{}N*".format(str(use_cycle)),basemask)
          
    # output folder
    outdir = os.path.join(rundir,"Cycle{}".format(str(use_cycle)))
    
    cmd = [bclToFastq, 
           "-b", 
           basemask,
           "-s",
           samplesheet,
           "-o",
           outdir,
           "-j",
           cores,
           "-t",
           "Cycle{}".format(str(use_cycle)),
           "-x"]
      
    # Call the base conversion
    output.append(subprocess.check_output(cmd))
    output.append("{}\tFinished processing run {} at cycle {}".format(datetime.datetime.utcnow().isoformat(),rundir,use_cycle))
      
    # Upload the results to the destination
    if destination:
        output.append("{}\tUploading to {}".format(datetime.datetime.utcnow().isoformat(),destination))
        cmd = ["rsync", 
               "--checksum", 
               "--archive",  
               "--prune-empty-dirs",
               "--include=*/",
               "--include={}**/*.xml".format(os.path.basename(outdir)),
               "--include={}**/*.htm".format(os.path.basename(outdir)),
               "--include={}**/*.csv".format(os.path.basename(outdir)),
               "--include={}**/*.fastq.gz".format(os.path.basename(outdir)),
               "--exclude=*",
               os.path.dirname(outdir), 
               destination
              ]
        output.append(subprocess.check_output(cmd))
        
        output.append("{}\tSetting permissions on {}".format(datetime.datetime.utcnow().isoformat(),destination))
        cmd = ["ssh", 
               destination.split(":")[0], 
               "chgrp",
               "-R",
               "uppmax",
               os.path.join(destination.split(":")[1],os.path.basename(rundir))
               ]  
        output.append(subprocess.check_output(cmd))
        cmd = ["ssh", 
               destination.split(":")[0], 
               "chmod",
               "-R",
               "g+rwX",
               os.path.join(destination.split(":")[1],os.path.basename(rundir))
               ]    
        output.append(subprocess.check_output(cmd))
        
        output.append("{}\tFinished uploading to {}".format(str(datetime.datetime.utcnow().isoformat()),destination))
        
    duration = (datetime.datetime.utcnow()-starttime).total_seconds()
    h = int(duration/3600)
    duration = duration%3600
    m = int(duration/60)
    duration = int(duration%60)
    output.append("{}\tAll finished at cycle {}, total duration was {}h, {}m, {}s".format(datetime.datetime.utcnow().isoformat(),
                                                                        use_cycle,
                                                                        str(h),
                                                                        str(m),
                                                                        str(duration)))

    # Send an email notification
    try:
        msg = MIMEText("\n".join(output))
        msg['To'] = ",".join(recipients)
        msg['Subject'] = "Cycle {} of run {} extracted and uploaded".format(use_cycle,os.path.basename(rundir))
        msg['From'] = sender
        s = smtplib.SMTP("localhost")
        s.sendmail(msg['From'], recipients, msg.as_string())
        s.quit()   
    except:
        print("\n".join(output))
        
