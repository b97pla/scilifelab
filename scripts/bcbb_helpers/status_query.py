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
            print "{}***ERROR***: Could not locate samplesheet in flowcell directory. Skipping.."
            continue
        fc_status['samplesheet'] = samplesheet

        # Get a list of the projects in the samplesheet
        projects = _get_projects(samplesheet,project)
        if len(projects) == 0:
            print "\t***WARNING***: No projects matched your filter [{}] for flowcell. Skipping..".format(project)
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
        print "{}".format(fc_status['flowcell'])
        for proj_status in fc_status.get('projects',[]):
            
            if brief:
                print "\t".join([proj_status['project'],
                                "All Finished" if proj_status.get('finished',False) else "{}/{} finished".format(proj_status.get('no_finished_samples',0),
                                                                                                                 len(proj_status.get('samples',[])))])
                continue
            
            print "\t{}".format(proj_status['project'])
            print "\t\tAnalysis finished...\t{}".format("Yes" if proj_status.get('finished',False) else "{}/{} samples".format(proj_status.get('no_finished_samples',0),
                                                                                                                               len(proj_status.get('samples',[]))))
            print "\t\tProject analysis folder exists...\t{}".format("Yes" if proj_status.get('project_dir',None) is not None else "No")
            
            for sample_status in proj_status.get('samples',[]):
                
                print "\t\t{}".format(sample_status['sample_id'])
                
                print "\t\t\tAnalysis finished...\t{}".format("Yes" if sample_status.get('finished',False) else "No") 
                stat, tstamp = sample_status.get('pipeline_progress',[None,None])
                print "\t\t\tMost recent pipeline step completed...\t{}".format("{}\t{}".format(os.path.splitext(os.path.basename(stat))[0],
                                                                                                tstamp.isoformat()) if tstamp is not None else "N/A")
                stat, tstamp = sample_status.get('pipeline_started',[None,None])
                print "\t\t\tPipeline started...\t{}".format("{}\t{}".format("Yes",
                                                                             tstamp.isoformat()) if tstamp is not None else "No")
                stat, tstamp = sample_status.get('pipeline_log',[None,None])
                print "\t\t\tPipeline log file last updated...\t{}".format(tstamp.isoformat() if tstamp is not None else "N/A")
                
                slurm_jobs = sample_status.get('slurm_job',[])
                print "\t\t\tPipeline jobs exist in slurm...\t{}".format("Yes" if len(slurm_jobs) > 0 else "No")
                for job in slurm_jobs:
                    print "\t\t\t\tJob id {} is {}".format(str(job[0]), str(job[1]))
                
                stat, finished = sample_status.get('fastq_screen',[None,False])
                print "\t\t\tFastq_screen has finished...\t{}".format("Yes" if finished else "{}".format("No" if stat is not None else "N/A"))
                print "\t\t\tSample flowcell analysis folder exists...\t{}".format("Yes" if sample_status.get('sample_fc_dir',None) is not None else "No")
                print "\t\t\tSample analysis folder exists...\t{}".format("Yes" if sample_status.get('sample_dir',None) is not None else "No")
                
                
            
        
        
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

### ---- Tests ----

import unittest
import tempfile
import shutil
import random
import string
import tests.generate_test_data as td
from mock import Mock

class TestMetadata(unittest.TestCase):
    
    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="status_query_metadata_test_")
        
    def tearDown(self):
        shutil.rmtree(self.rootdir)
        
        
    def test__get_flowcelldirs(self):
        """Test that the _get_flowcelldirs method behaves as expected
        """
        # Create a few random files and folders and assert that they are not returned
        for n in range(3):
            os.mkdir(os.path.join(self.rootdir,''.join(random.choice(string.ascii_uppercase) for x in range(5))))
            fh, _ = tempfile.mkstemp(dir=self.rootdir)
            os.close(fh)
            
        self.assertListEqual([], _get_flowcelldirs(self.rootdir),
                             "Listing available flowcells did not return an empty list")
        
        # Create a few flowcell directories and make sure that they are returned
        fcdirs = []
        for n in range(3):
            fcdirs.append(os.path.join(self.rootdir,td.generate_run_id()))
            os.mkdir(fcdirs[-1])
        
        fcdirs = sorted(fcdirs)
        self.assertListEqual(fcdirs,sorted(_get_flowcelldirs(self.rootdir)),
                             "Listing available flowcells did not return the 3 expected folders")
        
        # Test that filtering based on flowcell name works
        self.assertListEqual([fcdirs[-1]],_get_flowcelldirs(self.rootdir,os.path.basename(fcdirs[-1])),
                             "Filtering based on flowcell name did not returned the expected output")
        
        # Test that filtering based on partial flowcell name works
        self.assertListEqual([fcdirs[-1]],_get_flowcelldirs(self.rootdir,os.path.basename(fcdirs[-1][10:])),
                             "Filtering based on partial flowcell name did not returned the expected output")
    
    def test__get_samplesheet(self):
        """Test that the _get_samplesheet method behaves as expected
        """
        # Create a few random files and folders and assert that they are not returned
        suffixes = [".csv","",""]
        for n in range(3):
            os.mkdir(os.path.join(self.rootdir,''.join(random.choice(string.ascii_uppercase) for x in range(5))))
            fh, _ = tempfile.mkstemp(dir=self.rootdir, suffix=suffixes[n])
            os.close(fh)
            
        self.assertIsNone(_get_samplesheet(self.rootdir),
                          "Getting non-existing samplesheet did not return None")
        
        # Create a SampleSheet.csv and a [FCID].csv file and assert that they are
        # returned with a preference for the [FCID].csv file
        fcid = td.generate_fc_barcode()
        fcdir = os.path.join(self.rootdir,td.generate_run_id(fc_barcode=fcid))
        os.mkdir(fcdir)
        
        ss = [os.path.join(fcdir,"SampleSheet.csv"),
              os.path.join(fcdir,"{}.csv".format(fcid))]
        for s in ss:
            utils.touch_file(s)
            self.assertEqual(s,_get_samplesheet(fcdir),
                             "Did not get existing {}".format(os.path.basename(s)))
            
    def test__get_projects(self):
        """Test that getting the projects from a flowcell works
        """
        # Assert that an empty file returns an empty list
        fh, ssheet = tempfile.mkstemp(dir=self.rootdir, suffix=".csv")
        os.close(fh)
        self.assertListEqual([],_get_projects(ssheet),
                             "The list of projects for an empty file is not empty")
        
        # Generate artificial samplesheet data
        data = td.generate_samplesheet_data()
        projects = {}
        for d in data:
            projects[d[-1]] = 1
        
        # Write the data to a samplesheet
        td._write_samplesheet(data,ssheet)
         
        # Assert that the list of projects returned is the same that we generated
        self.assertListEqual(sorted(projects.keys()),sorted(_get_projects(ssheet)),
                             "The list of projects does not match the original list")
        
        # Assert that the list of projects returned is filtered as specified
        self.assertListEqual([projects.keys()[-1]],_get_projects(ssheet,projects.keys()[-1]),
                             "The filtered list of projects does not match the expected")
        
    def test__get_project_analysis_dir(self):
        """Test that getting the project analysis folder behaves as expected
        """
        # Assert that none is returned when no folder exists
        proj = td.generate_project()
        self.assertIsNone(_get_project_analysis_dir(self.rootdir,proj),
                          "Did not return empty result for non-existing folders")
        
        # Assert that none is still returned when some mismatching folders exist
        for n in range(5):
            os.mkdir(os.path.join(self.rootdir,td.generate_project()))
        self.assertIsNone(_get_project_analysis_dir(self.rootdir,proj),
                          "Did not return empty result for mismatching folders")
        
        # Assert that a file with the same name as the project is not returned
        projdir = os.path.join(self.rootdir,proj)
        utils.touch_file(projdir)
        self.assertIsNone(_get_project_analysis_dir(self.rootdir,proj),
                          "Returned a file with matching name. Should only return folders")
        os.unlink(projdir)
        
        # Assert that the corrct folder is returned when it exists
        os.mkdir(projdir)
        self.assertEqual(projdir,_get_project_analysis_dir(self.rootdir,proj),
                         "The expected project folder was not returned")
        
    def test__get_project_samples(self):
        """Test that getting the project samples from a samplesheet behaves as expected
        """
        
        # Generate artificial samplesheet data
        data = td.generate_samplesheet_data()
        fh, ssheet = tempfile.mkstemp(dir=self.rootdir, suffix=".csv")
        os.close(fh)
        td._write_samplesheet(data,ssheet)
         
        # Assert that getting samples for a non-existing project returns an empty list
        self.assertListEqual([],_get_project_samples(ssheet,td.generate_project()),
                             "Getting samples for a non-existing project returned unexpected output")
        
        # Iterate over the projects and assert that the returned samples are correct
        samples = {}
        for row in data:
            if row[9] not in samples:
                samples[row[9]] = []
            samples[row[9]].append(row[2])
        
        for proj, sample in samples.items():
            self.assertListEqual(sorted(sample),sorted(_get_project_samples(ssheet,proj)),
                                 "The returned list of samples did not match the original")

    def test__get_sample_log(self):
        """Test that getting the sample pipeline logfile behaves as expected
        """
        sample = td.generate_sample()
        # Assert that an empty directory returns None
        self.assertIsNone(_get_sample_log(self.rootdir,sample),
                          "Getting a sample log from an empty directory did not return None")
        
        # Assert that non-relevant log files are not returned
        for n in range(5):
            utils.touch_file(os.path.join(self.rootdir,"{}-bcbb.log".format(td.generate_sample())))
        self.assertIsNone(_get_sample_log(self.rootdir,sample),
                          "Getting a non-existing sample log did not return None")
        
        # Assert that the correct log file is returned when it exists
        utils.touch_file(os.path.join(self.rootdir,"{}-bcbb.log".format(sample)))
        self.assertEqual(os.path.join(self.rootdir,"{}-bcbb.log".format(sample)),_get_sample_log(self.rootdir,sample),
                         "Getting an existing sample log file did not return the expected file")
        
    def test__get_pipeline_indicator(self):
        """Getting pipeline indicator files
        """
        
        self.assertListEqual([],_get_pipeline_indicator(self.rootdir),
                          "Empty directory did not return an empty list")
        
        # Create some random files and assert that they are not picked up
        for n in range(5):
            fh, _ = tempfile.mkstemp(dir=self.rootdir)
            os.close(fh)
        self.assertListEqual([],_get_pipeline_indicator(self.rootdir),
                          "Non-existing indicator files did not return an empty list")
        
        # Create some indicator files and assert that they are returned
        ifiles = []
        for n in range(1,6):
            ifiles.append(os.path.join(self.rootdir,"{s:02d}_{act}.txt".format(s=n,act=''.join(random.choice(string.ascii_lowercase) for x in range(5)))))
            utils.touch_file(ifiles[-1])
            
        self.assertListEqual(sorted(ifiles),sorted(_get_pipeline_indicator(self.rootdir)),
                          "Existing indicator files did not return the expected output")
        
        # Assert that asking for a specific indicator returns the expected output
        self.assertListEqual(sorted(ifiles[0:2]),sorted(_get_pipeline_indicator(self.rootdir,range(1,3))),
                          "Specific indicator files did not return the expected output")
        
        # Assert that asking for a specific non-existing indicator returns the expected output
        self.assertListEqual([ifiles[-1]],sorted(_get_pipeline_indicator(self.rootdir,range(len(ifiles),len(ifiles)+2))),
                          "Specific non-existing indicator files did not return the expected output")
        
    def test__get_most_recent(self):
        """Get the most recent timestamp and corresponding file
        """
        
        # Assert that an empty input returns a (0,None) tuple
        self.assertEqual((datetime.datetime.fromtimestamp(0.0),None),_get_most_recent([]),
                             "No input files did not return expected output")
        
        # Assert that a file without timestamp returns (0,None)
        fh, tfile = tempfile.mkstemp(dir=self.rootdir)
        os.close(fh)
        self.assertEqual((datetime.datetime.fromtimestamp(0.0),None),_get_most_recent([tfile]),
                             "Empty input file did not return expected output")
        
        # Assert that the file with the most recent timestamp is returned
        tfiles = []
        for t in [1000.,2000.,3000.,4000.]:
             fh, tfile = tempfile.mkstemp(dir=self.rootdir)
             os.write(fh,"{}\n".format(datetime.datetime.fromtimestamp(t).isoformat()))
             os.close(fh)
             tfiles.append(tfile)
        
        self.assertEqual((datetime.datetime.fromtimestamp(4000.0),tfiles[-1]),_get_most_recent(tfiles),
                             "Timestamped input files did not return expected output")
        
        with open(tfiles[0],"a") as fh:
            fh.write("{}\n".format(datetime.datetime.fromtimestamp(5000.0).isoformat()))
        self.assertEqual((datetime.datetime.fromtimestamp(5000.0),tfiles[0]),_get_most_recent(tfiles),
                             "Input files with multiple timestamps did not return expected output")
            
    def test__fastq_screen_finished(self):
        """Detecting finished state of fastq_screen
        """
        
        # Assert that an empty directory doesn't indicate finished state
        self.assertFalse(_fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished without output files")
        
        # Create an output file and corresponding png but no rows in output
        sample_file = os.path.join(self.rootdir,"{}_fastq_screen.txt".format(td.generate_sample()))
        png_file = "{}.png".format(os.path.splitext(sample_file)[0])
        utils.touch_file(sample_file)
        utils.touch_file(png_file)
        self.assertFalse(_fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished with empty output file")
        
        # Write some output and assert fastq_screen is detected as finished
        with open(sample_file,"w") as fh:
            for n in range(5):
                fh.write("{}\n".format(str(n)))
        self.assertTrue(_fastq_screen_finished(self.rootdir),
                         "Fastq screen should be considered finished with non-empty output file and corresponding png")
        
        # Remove the png and assert fastq_screen is not finished
        os.unlink(png_file)
        self.assertFalse(_fastq_screen_finished(self.rootdir),
                         "Fastq screen should not be considered finished with non-empty output file but without corresponding png")
        
    def test__get_slurm_jobid(self):
        """Test that the jobid for a slurm job name can be extracted
        """
        
        # Mock the system calls
        subprocess.check_output = Mock(return_value='')
        # Assert that getting non-existing jobs return an empty job list
        self.assertListEqual([],_get_slurm_jobid("jobname"),
                             "Querying for jobid of non-existing job should return an empty list")
        # Assert that a returned job id is parsed correctly
        for jobids in [[123456789],[123456789,987654321]]:
            subprocess.check_output = Mock(return_value="\n".join([str(jid) for jid in jobids]))
            self.assertListEqual(jobids,_get_slurm_jobid("jobname"),
                                 "Querying for jobid of existing job did not return the correct value")
        
        
        
        
        
            
        
            
        