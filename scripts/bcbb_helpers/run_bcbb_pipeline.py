#!/usr/bin/env python
import drmaa
import os
import sys
import glob
import time
import yaml
import subprocess
import copy
import tempfile
import argparse

import bcbio.solexa.flowcell
import bcbio.solexa.samplesheet
from bcbio.pipeline.config_loader import load_config
import scilifelab.scripts.bcbb_helpers.report_to_gdocs as report

# The directory where CASAVA has written the demuxed output
CASAVA_OUTPUT_DIR = "Unaligned"
# The analysis script for running the pipeline in parallell mode (on one node)  
PARALLELL_ANALYSIS_SCRIPT="automated_initial_analysis.py"
# The analysis script for running the pipeline in distributed mode (across multiple nodes/cores)
DISTRIBUTED_ANALYSIS_SCRIPT="distributed_nextgen_pipeline.py"
# For non-CASAVA analysis, this script is used to sanitize the run_info.yaml configuration file
PROCESS_YAML_SCRIPT = "process_run_info.py"
# If True, will sanitize the run_info.yaml configuration file when running non-CASAVA analysis
PROCESS_YAML = True
# If True, will assign the distributed master process and workers to a separate RabbitMQ queue for each flowcell 
FC_SPECIFIC_AMPQ = True

def main(post_process_config_file, fc_dir, run_info_file=None, only_run=False, only_setup=False, ignore_casava=False):
    
    run_arguments = [[os.getcwd(),post_process_config_file,fc_dir,run_info_file]]
    if has_casava_output(fc_dir) and not ignore_casava:
        if not only_run:
            run_arguments = setup_analysis_directory_structure(post_process_config_file, fc_dir, run_info_file)
             
    else:
        if not only_run:
            run_arguments = setup_analysis(post_process_config_file, fc_dir, run_info_file)
    
    if not only_setup:
        for arguments in run_arguments:
            run_analysis(*arguments)
 
def run_analysis(work_dir, post_process, fc_dir, run_info):
    """Changes into the supplied work_dir directory and submits 
        the job using the supplied arguments and with slurm parameters
        obtained from the post_process.yaml configuration
    """
    
    # Move to the working directory
    start_dir = os.getcwd()
    os.chdir(work_dir)
    
    config = load_config(post_process)
    
    if str(config["algorithm"]["num_cores"]) == "messaging":
        analysis_script = DISTRIBUTED_ANALYSIS_SCRIPT
    else:
        analysis_script = PARALLELL_ANALYSIS_SCRIPT
        
    job_cl = [analysis_script, post_process, fc_dir, run_info]
    
    cp = config["distributed"]["cluster_platform"]
    cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
    platform_args = config["distributed"]["platform_args"].split()
    
    print "Submitting job"
    jobid = cluster.submit_job(platform_args, job_cl)
    print 'Your job has been submitted with id ' + jobid

    # Change back to the starting directory
    os.chdir(start_dir)

def setup_analysis(post_process_config, archive_dir, run_info_file):
    """Does a non-casava pre-analysis setup and returns a list of arguments
       that can be passed to the run_analysis function in order to start the
       analysis.
    """
    
    # Set the barcode type in run_info.yaml to "illumina", strip the 7th nucleotide and set analysis to 'Minimal'
    if run_info_file is not None and PROCESS_YAML:
        print "---------\nProcessing run_info:"
        run_info_backup = "%s.orig" % run_info_file
        os.rename(run_info_file,run_info_backup)
        cl = ["%s" % PROCESS_YAML_SCRIPT,run_info_backup,"--analysis","Align_illumina","--out_file",run_info_file,"--ascii","--clear_description"]
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
    base_dir = analysis["base_dir"]
    
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
        
        with open(post_process_config,"w") as fh:
            fh.write(yaml.safe_dump(config, default_flow_style=False, allow_unicode=True, width=1000)) 
            
    return [[os.getcwd(),post_process_config,archive_dir,run_info_file]]
        
def setup_analysis_directory_structure(post_process_config_file, fc_dir, custom_config_file):
    """Parse the CASAVA 1.8+ generated flowcell directory and create a 
       corresponding directory structure suitable for bcbb analysis,
       complete with sample-specific and project-specific configuration files.
       Returns a list of arguments, both sample- and project-specific, that can 
       be passed to the run_analysis method for execution
    """
    config = load_config(post_process_config_file)
    analysis_dir = os.path.abspath(config["analysis"]["base_dir"])
    assert os.path.exists(fc_dir), "ERROR: Flowcell directory %s does not exist" % fc_dir
    assert os.path.exists(analysis_dir), "ERROR: Analysis top directory %s does not exist" % analysis_dir
    
    # A list with the arguments to each run, when running by sample
    sample_run_arguments = []
    
    # Parse the flowcell dir
    fc_dir_structure = parse_casava_directory(fc_dir)
    [fc_date, fc_name] = [fc_dir_structure['fc_date'],fc_dir_structure['fc_name']]
    fc_run_id = "%s_%s" % (fc_date,fc_name)
    
    # Copy the basecall stats directory 
    _copy_basecall_stats(os.path.join(fc_dir_structure['fc_dir'],fc_dir_structure['basecall_stats_dir']), analysis_dir)
    
    # Parse the custom_config_file
    custom_config = []
    if custom_config_file is not None:
        with open(custom_config_file) as fh:
            custom_config = yaml.load(fh)
    
    # Iterate over the projects in the flowcell directory
    for project in fc_dir_structure.get('projects',[]):
        # Create a project directory if it doesn't already exist
        project_name = project['project_name'].replace('__','.')
        project_dir = os.path.join(analysis_dir,project_name)
        if not os.path.exists(project_dir):
            os.mkdir(project_dir,0770)
        
        src_project_dir = os.path.join(fc_dir_structure['fc_dir'],fc_dir_structure['data_dir'],project['project_dir'])
        project_samplesheet = _merge_samplesheets(src_project_dir, project_dir, project.get('samples',[]))
        
        # Iterate over the samples in the project
        for sample_no, sample in enumerate(project.get('samples',[])):
            # Create a directory for the sample if it doesn't already exist
            sample_name = sample['sample_name'].replace('__','.')
            sample_dir = os.path.join(project_dir,sample_name)
            if not os.path.exists(sample_dir):
                os.mkdir(sample_dir,0770)
            
            # Create a directory for the flowcell if it does not exist
            dst_sample_dir = os.path.join(sample_dir,fc_run_id)
            if not os.path.exists(dst_sample_dir):
                os.mkdir(dst_sample_dir,0770)
            
            # rsync the source files to the sample directory
            src_sample_dir = os.path.join(src_project_dir,sample['sample_dir'])
            sample_files = do_rsync([os.path.join(src_sample_dir,f) for f in (sample.get('files',[]) + sample.get('samplesheet',[]))],dst_sample_dir)
            
            # Generate a sample-specific configuration yaml structure
            samplesheet = os.path.join(src_sample_dir,sample['samplesheet'])
            sample_config = bcbb_configuration_from_samplesheet(samplesheet)
             
            # Append the sequence files to the config
            for lane in sample_config:
                if 'multiplex' in lane:
                    for sample in lane['multiplex']:
                        sample['files'] = [os.path.basename(f) for f in sample_files if f.find("_%s_L00%d_" % (sample['sequence'],int(lane['lane']))) >= 0]
                else:
                    lane['files'] = [os.path.basename(f) for f in sample_files if f.find("_L00%d_" % int(lane['lane'])) >= 0]
                    
            sample_config = override_with_custom_config(sample_config,custom_config)
            
            arguments = _setup_config_files(dst_sample_dir,sample_config,post_process_config_file,fc_dir_structure['fc_dir'],sample_name,fc_date,fc_name)
            sample_run_arguments.append([arguments[1],arguments[0],arguments[1],arguments[3]])
        
    return sample_run_arguments

def _merge_samplesheets(src_project_dir, dst_project_dir, samples):
    """Merge the samplesheets for a project's samples into one
    """
    # Parse the samplesheet content into a data structure
    ssheet_data = {}
    for sample in samples:
        sample_ssheet = sample.get('samplesheet',None)
        if sample_ssheet:
            with open(sample_ssheet) as in_handle:
                for row in in_handle:
                    ssheet_data[row.split(",")] = 1
    samplesheet = os.path.join(dst_project_dir, "SampleSheet.csv")
    with open(samplesheet,"w") as out_handle:
        for row in sorted(ssheet_data.keys(), key=lambda x: x[1]):
            out_handle.write(",".join(row))
            out_handle.write("\n")
    
    return samplesheet
                    
    
def _copy_basecall_stats(source_dir, destination_dir):
    """Copy relevant files from the Basecall_Stats_FCID directory
       to the analysis directory
    """
    
    # First create the directory in the destination
    dirname = os.path.join(destination_dir,os.path.basename(source_dir))
    try:
        os.mkdir(dirname)
    except:
        pass
    
    # List the files/directories to copy
    files = glob.glob(os.path.join(source_dir,"*.htm"))
    files += glob.glob(os.path.join(source_dir,"*.xml"))
    files += glob.glob(os.path.join(source_dir,"*.xsl"))
    files += [os.path.join(source_dir,"Plots")]
    files += [os.path.join(source_dir,"css")]
    do_rsync(files,dirname)
 
def override_with_custom_config(org_config, custom_config):
    """Override the default configuration from the .csv samplesheets
       with a custom configuration. Will replace overlapping options
       or add options that are missing from the samplesheet-generated
       config.
    """
    
    new_config = copy.deepcopy(org_config)
    
    for item in new_config:
        for custom_item in custom_config:
            if item['lane'] != custom_item.get('lane',""):
                continue
            for key, val in custom_item.items():
                if key == 'multiplex':
                    continue
                item[key] = val
                
            for sample in item.get('multiplex',[]):
                if 'sequence' not in sample:
                    continue
                for custom_sample in custom_item.get('multiplex',[]):
                    if sample['sequence'] == custom_sample.get('sequence',""):
                        for key, val in custom_sample.items():
                            sample[key] = val
                        break
            break
        
    return new_config
       
def _setup_config_files(dst_dir,configs,post_process_config_file,fc_dir,sample_name="run",fc_date=None,fc_name=None):
    
    # Setup the data structure
    config_data_structure = {'details': configs}
    if fc_date is not None:
        config_data_structure['fc_date'] = fc_date
    if fc_name is not None:
        config_data_structure['fc_name'] = fc_name
        
    # Dump the config to file
    config_file = os.path.join(dst_dir,"%s-bcbb-config.yaml" % sample_name)
    with open(config_file,'w') as fh:
        fh.write(yaml.safe_dump(config_data_structure, default_flow_style=False, allow_unicode=True, width=1000))
            
    # Copy post-process file
    with open(post_process_config_file) as fh:
        local_post_process = yaml.load(fh) 
    # Update galaxy config to point to the original location
    local_post_process['galaxy_config'] = bcbio.utils.add_full_path(local_post_process['galaxy_config'],os.path.abspath(os.path.dirname(post_process_config_file)))
    # Add job name and output paths to the cluster platform arguments
    if 'distributed' in local_post_process and 'platform_args' in local_post_process['distributed']:
        slurm_out = "%s-bcbb.log" % sample_name
        local_post_process['distributed']['platform_args'] = "%s -J %s -o %s -D %s" % (local_post_process['distributed']['platform_args'], sample_name, slurm_out, dst_dir)
            
    local_post_process_file = os.path.join(dst_dir,"%s-post_process.yaml" % sample_name)
    with open(local_post_process_file,'w') as fh:
        fh.write(yaml.safe_dump(local_post_process, default_flow_style=False, allow_unicode=True, width=1000))
            
    # Write the command for running the pipeline with the configuration files
    run_command_file = os.path.join(dst_dir,"%s-bcbb-command.txt" % sample_name)
    with open(run_command_file,"w") as fh:
        fh.write(" ".join([os.path.basename(__file__),"--only-run",os.path.basename(local_post_process_file), os.path.join("..",os.path.basename(dst_dir)), os.path.basename(config_file)])) 
        fh.write("\n")   
    
    return [os.path.basename(local_post_process_file), dst_dir, fc_dir, os.path.basename(config_file)]
    
def bcbb_configuration_from_samplesheet(csv_samplesheet):
    """Parse an illumina csv-samplesheet and return a dictionary suitable for the bcbb-pipeline
    """
    tfh, yaml_file = tempfile.mkstemp('.yaml','samplesheet')
    os.close(tfh)
    yaml_file = bcbio.solexa.samplesheet.csv2yaml(csv_samplesheet,yaml_file)
    with open(yaml_file) as fh:
        config = yaml.load(fh)
    
    # Replace the default analysis
    ## TODO: This is an ugly hack, should be replaced by a custom config 
    for lane in config:
        if lane.get('genome_build','') == 'hg19':
            lane['analysis'] = 'Align_standard_seqcap'
        else:
            lane['analysis'] = 'Align_standard'
        for plex in lane.get('multiplex',[]):
            if plex.get('genome_build','') == 'hg19':
                plex['analysis'] = 'Align_standard_seqcap'
            else:
                plex['analysis'] = 'Align_standard'
                
    # Remove the yaml file, we will write a new one later
    os.remove(yaml_file)
    
    return config
                
def do_rsync(src_files, dst_dir):
    cl = ["rsync","-car"]
    cl.extend(src_files)
    cl.append(dst_dir)
    cl = [str(i) for i in cl]
    # For testing, just touch the files rather than copy them
    # for f in src_files:
    #    open(os.path.join(dst_dir,os.path.basename(f)),"w").close()
    subprocess.check_call(cl)
    
    return [os.path.join(dst_dir,os.path.basename(f)) for f in src_files]
        
def parse_casava_directory(fc_dir):
    """Traverse a CASAVA 1.8+ generated directory structure and return a dictionary
    """ 
    projects = []
    
    fc_dir = os.path.abspath(fc_dir)
    fc_name, fc_date = bcbio.solexa.flowcell.get_flowcell_info(fc_dir)
    unaligned_dir = os.path.join(fc_dir,CASAVA_OUTPUT_DIR)
    basecall_stats_dir_pattern = os.path.join(unaligned_dir,"Basecall_Stats_*")
    basecall_stats_dir = None
    try:
        basecall_stats_dir = os.path.relpath(glob.glob(basecall_stats_dir_pattern)[0],fc_dir)
    except:
        print "WARNING: Could not locate basecall stats directory under %s" % unaligned_dir
        
    project_dir_pattern = os.path.join(unaligned_dir,"Project_*")
    for project_dir in glob.glob(project_dir_pattern):
        project_samples = []
        sample_dir_pattern = os.path.join(project_dir,"Sample_*")
        for sample_dir in glob.glob(sample_dir_pattern):
            fastq_file_pattern = os.path.join(sample_dir,"*.fastq.gz")
            samplesheet_pattern = os.path.join(sample_dir,"*.csv")
            fastq_files = [os.path.basename(file) for file in glob.glob(fastq_file_pattern)]
            samplesheet = glob.glob(samplesheet_pattern)
            assert len(samplesheet) == 1, "ERROR: Could not unambiguously locate samplesheet in %s" % sample_dir
            sample_name = sample_dir.replace(sample_dir_pattern[0:-1],'')
            project_samples.append({'sample_dir': os.path.relpath(sample_dir,project_dir), 'sample_name': sample_name, 'files': fastq_files, 'samplesheet': os.path.basename(samplesheet[0])})
        project_name = project_dir.replace(project_dir_pattern[0:-1],'')
        projects.append({'project_dir': os.path.relpath(project_dir,unaligned_dir), 'project_name': project_name, 'samples': project_samples})
    
    return {'fc_dir': fc_dir, 'fc_name': fc_name, 'fc_date': fc_date, 'data_dir': os.path.relpath(unaligned_dir,fc_dir), 'basecall_stats_dir': basecall_stats_dir, 'projects': projects}
    
def has_casava_output(fc_dir):
    try:
        structure = parse_casava_directory(fc_dir)
        if len(structure['projects']) > 0:
            return True
    except:
        pass
    return False

def report_to_gdocs(fc_dir, post_process_config_file):
    # Rename any existing run_info.yaml as it will interfere with gdocs upload
    run_info = os.path.join(fc_dir, "run_info.yaml")
    if os.path.exists(run_info):
        os.rename(run_info, "{}.bak".format(run_info))
    report.main(os.path.basename(os.path.abspath(fc_dir)), post_process_config_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Wrapper script for bcbb pipeline. If given a .yaml configuration file, "\
                                     "a run folder containing the sequence data from an illumina run and, optionally, "\
                                     "a custom yaml file with options that should override what is specified in the "\
                                     "config and samplesheet, the script will copy the relevant files from [store_dir] "\
                                     "to [base_dir] and submit the automated_initial_analysis.py pipeline script for each "\
                                     "sample to the cluster platform specified in the configuration.")

    parser.add_argument("config", action="store", default=None, help="Path to the .yaml pipeline configuration file")
    parser.add_argument("fcdir", action="store", default=None, help="Path to the archive run folder")
    parser.add_argument("custom_config", action="store", default=None, help="Path to a custom configuration file with lane or sample specific options that will override the main configuration", nargs="?")
    parser.add_argument("-r", "--only-run", dest="only_run", action="store_true", default=False, help="Don't setup the analysis directory, just start the pipeline")
    parser.add_argument("-s", "--only-setup", dest="only_setup", action="store_true", default=False, help="Setup the analysis directory but don't start the pipeline")
    parser.add_argument("-i", "--ignore-casava", dest="ignore_casava", action="store_true", default=False, help="Ignore any Casava 1.8+ file structure and just assume the pre-casava pipeline setup")
    parser.add_argument("-g", "--no-google-report", dest="no_google_report", action="store_true", default=False, help="Don't upload any demultiplex statistics to Google Docs")
    args = parser.parse_args()
        
    main(args.config,args.fcdir,args.custom_config,args.only_run,args.only_setup,args.ignore_casava)
    if not args.no_google_report:
        report_to_gdocs(args.fcdir, args.config)

# --- Testing code: run with 'nosetests -v -s run_bcbb_pipeline.py'

import unittest
import tempfile
import random
import date
import shutil
import bcbio.utils as utils

def generate_fc_barcode():
    """Generate a flowcell barcode on the format ABC123CXX
    """
    return "".join(random.choice("{}{}CXX".format(string.ascii_lowercase,string.digits)) for i in xrange(6))

def generate_run_id(fc_barcode=generate_fc_barcode()):
    """Generate a run identifier
    """
    return "{}_SN{}_0{}_{}{}".format(datetime.date.today().strftime("%y%m%d"),
                                           random.randint(101,9999),
                                           random.randint(101,999),
                                           random.choice("AB"),
                                           fc_barcode)

def generate_project():
    """Generate a project name on the format J.Doe_90_11
    """
    return "{}.{}_{}_{}".format(
        random.choice(string.ascii_uppercase),
        "".join(random.choice(string.ascii_lowercase) for i in xrange(6)).capitalize(),
        random.randint(90,99), 
        random.randint(11,99))

def generate_sample(project_id=random.randint(500,999)):
    """Generate a sample name on the format P123_123B
    """
    return "P{}_{}{}_index{}".format(
                    project_id,
                    random.randint(101,199),
                    ["","A","B","C","F"][random.randint(0,4)],
                    random.randint(1,24))

def generate_sample_file(sample_name=generate_sample(), barcode=generate_barcode(), lane=random.randint(1,8), readno=1):
    """Generate a Casava 1.8+-style sample file name
    """
    return "{}_{}_L00{}_R{}_001.fastq.gz".format(sample_name,
                                                barcode,
                                                lane,
                                                readno)

def generate_barcode(len=6):
    """Generate a nucleotide barcode
    """
    return (random.choice("ACGT") for i in xrange(len))

def generate_run_samplesheet(barcode=generate_bc_barcode(), dst_file=None):
    
    csv_data = []
    
    # Generate data for each lane
    for lane in range(1,9):
        # Create some projects
        for i in range(3):
            project_name = generate_project()
            # Create some samples
            for j in range(3):
                sample_name = generate_sample()
                # Append the generated data to the csv_data
                csv_data.append([barcode,
                                 lane,
                                 sample_name,
                                 'hg19',
                                 generate_barcode(),
                                 project_name,
                                 'C',
                                 'R',
                                 'O',
                                 project_name])
    return _write_samplesheet(csv_data,dst_file)

def _write_samplesheet(csv_data, dst_file=None):
    
    header = ["FCID",
              "Lane",
              "SampleID",
              "SampleRef",
              "Index",
              "Description",
              "Control",
              "Recipe",
              "Operator",
              "SampleProject"]
    
    if dst_file is None:
        dst_file = tempfile.mkstemp(suffix=".csv", prefix="SampleSheet")
    with open(dst_file,"w") as out_handle:
        for row in (header + sorted(csv_data, key=lambda x: x[1])):
            out_handle.write(",".join(row))
            out_handle.write("\n")
    return dst_file
        

class CasavaStructureTest(unittest.TestCase):
    """Test for the Casava file structure
    """
    
    @classmethod
    def setUpClass(self):
        
        # Create the file structure in a temporary location
        self.test_archive_dir = tempfile.mkdtemp()
        
        # Create a flowcell id
        barcode = generate_bc_barcode()
        fcid = generate_run_id(fc_barcode=barcode)
        self.fc_dir = os.path.join(self.test_archive_dir,fcid)
        os.makedirs(self.fc_dir)
        
        # Generate run data and write it to a samplesheet
        self.samplesheet = os.path.join(self.fc_dir,"SampleSheet.csv")
        generate_run_samplesheet(barcode,self.samplesheet)
        
        # Create the file structure according to the samplesheet
        with open(self.samplesheet) as in_handle:
            for row in in_handle:
                if len(row) == 0 or row[0] == "#":
                    continue
                csv_data = row.split(",")
                lane = 0
                try:
                    lane = int(str(csv_data[1]))
                except ValueError:
                    # This is most likely the header row
                    pass
                
                sample_name = csv_data[2]
                project_name = csv_data[9]
                index_sequence = csv_data[4]
                
                project_folder = os.path.join(self.fc_dir,"Unaligned","Project_{}".format(project_name))
                sample_folder = os.path.join(project_folder,"Sample_{}".format(sample_name))
                sample_file_r1 = os.path.join(sample_folder,"{}_{}_L00{}_R1_001.fastq.gz".format(
                                                    sample_name,
                                                    index_sequence,
                                                    lane))
                sample_file_r2 = sample_file_r1.replace("_R1_","_R2_")
                sample_ssheet = os.path.join(sample_folder,"SampleSheet.csv")
                
                os.makedirs(sample_folder)
                utils.touch_file(sample_file_r1)
                utils.touch_file(sample_file_r2)
                sample_ssheet = _write_samplesheet(csv_data,sample_ssheet)
                
        # Create a Basecall_Stats_[FCID] folder and a Demultiplex_Stats.htm file
        bcall_dir = os.path.join(self.fc_dir,"Unaligned","Basecall_Stats_{}".format(barcode))
        os.mkdir(bcall_dir)
        utils.touch_file(os.path.join(bcall_dir,"Demultiplex_Stats.htm"))
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.test_archive_dir)
        
    def test_has_casava_output(self):
        self.assertTrue(has_casava_output(self.fc_dir))
            
class RunInfoYamlTest(unittest.TestCase):
    
    def setUp(self):
        """Set up a samplesheet to process
        """
        self.samplesheet = tempfile.mkstemp(suffix=".csv", prefix="SampleSheet")
        content = []
        
        # Create a barcode
        barcode = generate_barcode()
        
        # Create some projects
        projects = [generate_project() for i in range(5)]
        # Distribute projects and samples among the first 6 lanes
        for lane in range(1,7):
            # Pick a number of projects for each lane
            for i in range(random.randint(0,len(projects))):
                # Generate a number of samples for each project
                for j in range(random.randint(1,5)):
                    row = [barcode,lane,generate_sample(),generate_barcode(),"","","",projects[i]]
                    content.append(row)
        # Append one non-demultiplexed lane with a project that has samples in a multiplexed lane
        content.append([barcode,7,generate_sample(),"","","","",projects[0]])
        # Append a non-demultiplexed lane with a project without any other samples
        content.append([barcode,8,generate_sample(),"","","","",generate_project()]) 
        
        # Write the content to the samplesheet
        with open(self.samplesheet,"w") as out_handle:
            for row in content:
                out_handle.write(",".join(row))
                out_handle.write("\n")
                
    def tearDown(self):
        os.unlink(self.samplesheet)
        
