import glob
import os
import collections
import argparse
import csv
import bcbio.utils as utils
import datetime
from  dateutil  import  parser
import getpass
import subprocess
import drmaa

def status_query(archive_dir, analysis_dir, flowcell, project, brief):
    """Get a status report of the progress of flowcells based on a snapshot of the file system
    """
    
    last_step = 14
    status = []
    # Process each flowcell in the archive directory
    for fcdir in _get_flowcelldirs(archive_dir,flowcell):
        fc_status = {}
        fc_status['flowcell'] = os.path.basename(fcdir)
        
        # Locate the samplesheet
        samplesheet = _get_samplesheet(fcdir)
        if samplesheet is None:
            print("{}***ERROR***: Could not locate samplesheet in flowcell directory. Skipping..")
            continue
        fc_status['samplesheet'] = samplesheet

        # Get a list of the projects in the samplesheet
        projects = _get_projects(samplesheet,project)
        if len(projects) == 0:
            print("\t***WARNING***: No projects matched your filter [{}] for flowcell. Skipping..".format(project))
            continue
        
        fc_status['projects'] = []
        
        # Iterate over the projects in the flowcell
        for proj in projects:
            proj = proj.replace("__",".")
            proj_status = {}
            proj_status['project'] = proj
            
            pdir = _get_project_analysis_dir(analysis_dir, proj)
            if not pdir:
                continue
            
            proj_status['project_dir'] = pdir
            proj_status['samples'] = []
            proj_status['no_finished_samples'] = 0
            samples = _get_project_samples(samplesheet, proj)
            for smpl in samples:
                smpl = smpl.replace("__",".")
                sample_status = {}
                proj_status['samples'].append(sample_status)
                sample_status['sample_id'] = smpl
                sdir = _get_sample_analysis_dir(pdir, smpl)
                if not sdir:
                    continue
                sample_status['sample_dir'] = sdir
                
                # Match the flowcell we're processing to the sample flowcell directories
                sample_fc = [d for d in _get_flowcelldirs(sdir) if d.split("_")[-1] == fcdir.split("_")[-1]]
                if len(sample_fc) == 0:
                    continue
                sample_fc = sample_fc[0]
                sample_status['sample_fc_dir'] = sample_fc
                
                fastq_screen = _get_fastq_screen_folder(sample_fc)
                if fastq_screen:
                    sample_status['fastq_screen'] = [fastq_screen,_fastq_screen_finished(fastq_screen)]
                
                now = datetime.datetime.now()
                pipeline_start_indicator = _get_pipeline_indicator(sample_fc,[1])
                if len(pipeline_start_indicator) == 0:
                    continue
                pipeline_start_indicator = pipeline_start_indicator[0]
                
                most_recent, _ = _get_most_recent([pipeline_start_indicator])
                sample_status['pipeline_started'] = [pipeline_start_indicator,most_recent]
                
                most_recent, ifile = _get_most_recent(_get_pipeline_indicator(sample_fc))
                sample_status['pipeline_progress'] = [ifile,most_recent]
                
                sample_log = _get_sample_log(sample_fc,smpl)
                if not sample_log:
                    continue
                st = os.stat(sample_log)
                sample_status['pipeline_log'] = [sample_log,datetime.datetime.fromtimestamp(st.st_mtime)]
                
                jobids = _get_slurm_jobid(smpl)
                sample_status['slurm_job'] = []
                for jobid in jobids:
                    sample_status['slurm_job'].append([jobid,_get_slurm_jobstatus(jobid)])
                
                most_recent, ifile = _get_most_recent(_get_pipeline_indicator(sample_fc,[last_step]))
                if ifile is not None and sample_status.get('fastq_screen',[None,False])[1]:
                    sample_status['finished'] = True
                    proj_status['no_finished_samples'] += 1
                
            
            if proj_status['no_finished_samples'] == len(samples):
                proj_status['finished'] = True
                
            fc_status['projects'].append(proj_status)
            
        status.append(fc_status) 
    print_status(status,brief)
         
def print_status(status, brief=False):
    """Pretty-print the status output
    """
    
    for fc_status in status:
        print("{}".format(fc_status['flowcell']))
        for proj_status in fc_status.get('projects',[]):
            
            if brief:
                print("\t".join([proj_status['project'],
                                "All Finished" if proj_status.get('finished',False) else "{}/{} finished".format(proj_status.get('no_finished_samples',0),
                                                                                                                 len(proj_status.get('samples',[])))]))
                continue
            
            print("\t{}".format(proj_status['project']))
            print("\t\tAnalysis finished...\t{}".format("Yes" if proj_status.get('finished',False) else "{}/{} samples".format(proj_status.get('no_finished_samples',0),
                                                                                                                               len(proj_status.get('samples',[])))))
            print("\t\tProject analysis folder exists...\t{}".format("Yes" if proj_status.get('project_dir',None) is not None else "No"))
            
            for sample_status in proj_status.get('samples',[]):
                
                print("\t\t{}".format(sample_status['sample_id']))
                
                print("\t\t\tAnalysis finished...\t{}".format("Yes" if sample_status.get('finished',False) else "No") )
                stat, tstamp = sample_status.get('pipeline_progress',[None,None])
                print("\t\t\tMost recent pipeline step completed...\t{}".format("{}\t{}".format(os.path.splitext(os.path.basename(stat))[0],
                                                                                                tstamp.isoformat()) if tstamp is not None else "N/A"))
                stat, tstamp = sample_status.get('pipeline_started',[None,None])
                print("\t\t\tPipeline started...\t{}".format("{}\t{}".format("Yes",
                                                                             tstamp.isoformat()) if tstamp is not None else "No"))
                stat, tstamp = sample_status.get('pipeline_log',[None,None])
                print("\t\t\tPipeline log file last updated...\t{}".format(tstamp.isoformat() if tstamp is not None else "N/A"))
                
                slurm_jobs = sample_status.get('slurm_job',[])
                print("\t\t\tPipeline jobs exist in slurm...\t{}".format("Yes" if len(slurm_jobs) > 0 else "No"))
                for job in slurm_jobs:
                    print("\t\t\t\tJob id {} is {}".format(str(job[0]), str(job[1])))
                
                stat, finished = sample_status.get('fastq_screen',[None,False])
                print("\t\t\tFastq_screen has finished...\t{}".format("Yes" if finished else "{}".format("No" if stat is not None else "N/A")))
                print("\t\t\tSample flowcell analysis folder exists...\t{}".format("Yes" if sample_status.get('sample_fc_dir',None) is not None else "No"))
                print("\t\t\tSample analysis folder exists...\t{}".format("Yes" if sample_status.get('sample_dir',None) is not None else "No"))
                
                
            
        
        
def _fastq_screen_finished(fastq_screen_dir):
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

def _get_fastq_screen_folder(analysis_dir):
    """Get the fastq_screen output folder
    """
    return _get_dir(analysis_dir,"fastq_screen")

def _get_flowcelldirs(archive_dir, flowcell_id="CXX"):
    """List the available flowcells in the archive dir, optionally filtering on flowcell_id
    """
    pattern = os.path.join(archive_dir,"*{}".format(flowcell_id))
    matches = []
    for d in glob.glob(pattern):
        if not os.path.isdir(d):
            continue
        matches.append(d)
        
    return matches

def _get_most_recent(ifiles):
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

def _get_project_analysis_dir(analysis_dir, project):
    """Return the analysis project dir if it exists, otherwise return none
    """
    return _get_dir(analysis_dir,project)

def _get_sample_analysis_dir(analysis_dir, sample):
    """Return the analysis sample dir if it exists, otherwise return none
    """
    return _get_dir(analysis_dir,sample)

def _get_sample_log(sample_dir, sample):
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

def _get_pipeline_indicator(sample_dir, steps=[]):
    """Get the pipeline indicator files
    """
    ifiles = []
    if len(steps) == 0:
        pattern = os.path.join(sample_dir,"[0-9][0-9]_*.txt")
        return [f for f in glob.glob(pattern) if os.path.isfile(f)]
    for step in steps:
        pattern = os.path.join(sample_dir,"{s:02d}_*.txt".format(s=step))
        ifiles += [f for f in glob.glob(pattern) if os.path.isfile(f)]
    return ifiles

def _get_project_samples(samplesheet, project):
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

def _get_projects(samplesheet, project=None):
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
    
def _get_samplesheet(flowcell_dir):
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

def _get_slurm_jobid(jobname,user=getpass.getuser()):
    """Attempt to get the job id for a slurm job name
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

def _get_slurm_jobstatus(jobid):
    """Get the status for a jobid
    """
    s = drmaa.Session()
    s.initialize()
    status = s.jobStatus(str(jobid))
    s.exit()
    return status
    

def main():
    
    parser = argparse.ArgumentParser(description="Query the status of flowcells, projects, samples that are organized "\
                                     "according to the CASAVA file structure")

    parser.add_argument('-f','--flowcell', action='store', default=None, 
                        help="only display the status for the specified flowcell")
    parser.add_argument('-p','--project', action='store', default=None, 
                        help="display only the status for the specified project")
    parser.add_argument('-b','--brief', action='store_true', default=False, 
                        help="display only a brief summary of the status status")
    parser.add_argument('-r','--archive-dir', dest='archive_dir', action='store', default="/proj/a2010002/archive", 
                        help="path to the folder containing flowcell data")
    parser.add_argument('-a','--analysis-dir', dest='analysis_dir', action='store', default="/proj/a2010002/nobackup/illumina", 
                        help="path to the folder containing project analysis data")
    
    args = parser.parse_args()
    status_query(args.archive_dir,args.analysis_dir,args.flowcell,args.project,args.brief)
      
if __name__ == "__main__":
    main()
        