
import os
import argparse
import datetime
import scilifelab.bcbio.filesystem as bcbio
import scilifelab.utils.slurm as slurm
from scilifelab.illumina import IlluminaRun
from scilifelab.illumina.hiseq import HiSeqRun

def status_query(archive_dir, analysis_dir, flowcell, project, brief):
    """Get a status report of the progress of flowcells based on a snapshot of the file system
    """
    
    last_step = 14
    status = []
    # Process each flowcell in the archive directory
    for fcdir in IlluminaRun.get_flowcell(archive_dir,flowcell):
        fc_status = {}
        fc_status['flowcell'] = os.path.basename(fcdir)
        
        # Locate the samplesheet
        samplesheet = IlluminaRun.get_samplesheet(fcdir)
        if samplesheet is None:
            print("{}***ERROR***: Could not locate samplesheet in flowcell directory. Skipping..")
            continue
        fc_status['samplesheet'] = samplesheet

        # Get a list of the projects in the samplesheet
        projects = HiSeqRun.get_project_names(samplesheet)
        if len(projects) == 0:
            print("\t***WARNING***: No projects matched your filter [{}] for flowcell. Skipping..".format(project))
            continue
        
        fc_status['projects'] = []
        
        # Iterate over the projects in the flowcell
        for proj in projects:
            proj = proj.replace("__",".")
            proj_status = {}
            proj_status['project'] = proj
            
            pdir = bcbio.get_project_analysis_dir(analysis_dir, proj)
            if not pdir:
                continue
            
            proj_status['project_dir'] = pdir
            proj_status['samples'] = []
            proj_status['no_finished_samples'] = 0
            samples = HiSeqRun.get_project_sample_ids(samplesheet, proj)
            for smpl in samples:
                smpl = smpl.replace("__",".")
                sample_status = {}
                proj_status['samples'].append(sample_status)
                sample_status['sample_id'] = smpl
                sdir = bcbio.get_sample_analysis_dir(pdir, smpl)
                if not sdir:
                    continue
                sample_status['sample_dir'] = sdir
                
                # Match the flowcell we're processing to the sample flowcell directories
                sample_fc = [d for d in IlluminaRun.get_flowcell(sdir) if d.split("_")[-1] == fcdir.split("_")[-1]]
                if len(sample_fc) == 0:
                    continue
                sample_fc = sample_fc[0]
                sample_status['sample_fc_dir'] = sample_fc
                
                fastq_screen = bcbio.get_fastq_screen_folder(sample_fc)
                if fastq_screen:
                    sample_status['fastq_screen'] = [fastq_screen,bcbio.fastq_screen_finished(fastq_screen)]
                
                now = datetime.datetime.now()
                pipeline_start_indicator = bcbio.get_pipeline_indicator(sample_fc,[1])
                if len(pipeline_start_indicator) == 0:
                    continue
                pipeline_start_indicator = pipeline_start_indicator[0]
                
                most_recent, _ = bcbio.get_most_recent_indicator([pipeline_start_indicator])
                sample_status['pipeline_started'] = [pipeline_start_indicator,most_recent]
                
                most_recent, ifile = bcbio.get_most_recent_indicator(bcbio.get_pipeline_indicator(sample_fc))
                sample_status['pipeline_progress'] = [ifile,most_recent]
                
                sample_log = bcbio.get_sample_pipeline_log(sample_fc,smpl)
                if not sample_log:
                    continue
                st = os.stat(sample_log)
                sample_status['pipeline_log'] = [sample_log,datetime.datetime.fromtimestamp(st.st_mtime)]
                
                jobids = slurm.get_slurm_jobid(smpl)
                sample_status['slurm_job'] = []
                for jobid in jobids:
                    sample_status['slurm_job'].append([jobid,slurm.get_slurm_jobstatus(jobid)])
                
                most_recent, ifile = bcbio.get_most_recent_indicator(bcbio.get_pipeline_indicator(sample_fc,[last_step]))
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
        