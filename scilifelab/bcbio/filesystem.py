"""Helper functions to locate and manipulate files in the bcbio analysis file structure
"""

import os
import glob
import datetime
import csv
from  dateutil  import  parser

def fastq_screen_finished(fastq_screen_dir):
    """Determine if the finished output from fastq_screen exists
    """
    pattern = os.path.join(fastq_screen_dir,"*_fastq_screen.txt")
    tfiles = glob.glob(pattern)
    if len(tfiles) == 0:
        return False
    
    for tf in tfiles:
        # Check that a corresponding png outfile exists
        if not os.path.exists("{}.png".format(os.path.splitext(tf)[0])):
            return False
        
        # Check that output exists beyond the header row
        rows = 0
        with open(tf) as fh:
            for row in fh:
                rows += 1
        
        if rows <= 1:
            return False
    
    return True

def get_fastq_screen_folder(analysis_dir):
    """Get the fastq_screen output folder
    """
    return _get_dir(analysis_dir,"fastq_screen")

def get_flowcelldirs(archive_dir, flowcell_id="CXX"):
    """List the available flowcells in the archive dir, optionally filtering on flowcell_id
    """
    pattern = os.path.join(archive_dir,"*{}".format(flowcell_id))
    matches = []
    for d in glob.glob(pattern):
        if not os.path.isdir(d):
            continue
        matches.append(d)
        
    return matches

def get_most_recent_indicator(ifiles):
    """Return a tuple with the most recent timestamp and the file that contains it
    """
    most_recent = (datetime.datetime.fromtimestamp(0.0),None)
    for ifile in ifiles:
        with open(ifile) as fh:
            for line in fh:
                try:
                    time = parser.parse(line.strip())
                    if time > most_recent[0]:
                        most_recent = (time,ifile)
                except ValueError:
                    # If we cannot parse time stamp, we may have grabbed some non-relevant file, in which case we'll break
                    break 
    return most_recent

def get_project_analysis_dir(analysis_dir, project):
    """Return the analysis project dir if it exists, otherwise return none
    """
    return _get_dir(analysis_dir,project)

def get_sample_analysis_dir(analysis_dir, sample):
    """Return the analysis sample dir if it exists, otherwise return none
    """
    return _get_dir(analysis_dir,sample)

def get_sample_pipeline_log(sample_dir, sample):
    """Return the log file where the pipeline writes output
    """
    logfile = os.path.join(sample_dir,"{}-bcbb.log".format(sample))
    if os.path.exists(logfile) and os.path.isfile(logfile):
        return logfile
    return None

def _get_dir(parent, dir):
    pdir = os.path.join(parent,dir)
    if os.path.exists(pdir) and os.path.isdir(pdir):
        return pdir
    return None

def get_pipeline_indicator(sample_dir, steps=[]):
    """Get the pipeline indicator files
    """
    ifiles = []
    if len(steps) == 0:
        pattern = os.path.join(sample_dir,"[0-9][0-9]_*.txt")
        ifiles = [f for f in glob.glob(pattern) if os.path.isfile(f)]
    for step in steps:
        pattern = os.path.join(sample_dir,"{s:02d}_*.txt".format(s=step))
        ifiles += [f for f in glob.glob(pattern) if os.path.isfile(f)]
    
    # Make sure that the located files are indicator files
    return [ifile for ifile in ifiles if _is_indicator(ifile)]

def _is_indicator(fname):
    """Parse a potential indicator file and return True if it only contains timestamps.
       Otherwise returns False
    """
    try:
        with open(fname) as fh:
            for line in fh:
                # Skip empty lines
                if len(line.strip()) == 0:
                    continue
                parser.parse(line.strip())
    except:
        return False
    return True
    
def get_project_samples(samplesheet, project):
    """Return the samples listed in the samplesheet for a project
    """
    samples = []
    with open(samplesheet) as fh:
        csvread = csv.DictReader(fh, dialect='excel')
        for row in csvread:
            p = row.get("SampleProject",None)
            if p is not None and (p == project or p.replace("__",".") == project):
                samples.append(row.get("SampleID"))
    return samples

def get_projects(samplesheet, project=None):
    """List the projects available in the samplesheet. Optionally filter by project name.
    """
    projects = {}
    with open(samplesheet) as fh:
        csvread = csv.DictReader(fh, dialect='excel')
        for row in csvread:
            p = row.get("SampleProject",None)
            if p is not None and (project is None or p == project or p.replace("__",".") == project):
                projects[p] = 1
    return sorted(projects.keys())
    
def get_samplesheet(flowcell_dir):
    """Get the samplesheet from the flowcell directory, returning firstly [FCID].csv and secondly SampleSheet.csv
    """
    pattern = os.path.join(flowcell_dir,"*.csv")
    ssheet = None
    for f in glob.glob(pattern):
        if not os.path.isfile(f):
            continue
        name, _ = os.path.splitext(os.path.basename(f))
        if flowcell_dir.endswith(name):
            return f
        if name == "SampleSheet":
            ssheet = f
    return ssheet
