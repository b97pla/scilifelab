# A script to help doing the deliveries.
# Now using the Casava directory structure
# The user is asked to provide a project ID, a run name, and an UPPMAX project

import sys
import os
import glob
import re
from datetime import datetime
import argparse
import stat
from subprocess import check_call, CalledProcessError
from scilifelab.utils.misc import filtered_walk, query_yes_no
from scilifelab.utils.timestamp import utc_time

def fixProjName(pname):
    newname = pname[0].upper()
    postperiod = False
    for i in range(1, len(pname)):
        if pname[i] == ".": 
            newname += pname[i]
            postperiod = True
        elif postperiod: 
            newname += pname[i].upper()
            postperiod = False
        else:
            newname += pname[i]
            postperiod = False
    return newname

def is_fastq(fname):
    fastq_ext = [".fastq.gz",
                 ".fastq",
                 "_fastq.txt.gz",
                 "_fastq.txt"]
    for ext in fastq_ext:
        if fname.endswith(ext):
            return True
    return False
        
def get_file_copy_list(proj_base_dir, dest_proj_path, fcid, deliver_all_fcs, deliver_nophix, skip_list, logfile, dry):
    to_copy = []
    for fqfile in filtered_walk(proj_base_dir, 
                                is_fastq, 
                                include_dirs=[fcid] if not deliver_all_fcs else None, 
                                exclude_dirs=skip_list):
        
        # Get the run_name and sample_name from the path
        sample_name, run_name, _ = os.path.relpath(fqfile,proj_base_dir).split(os.sep,2)
        date, fc_id = run_name.split('_')
            
        # Skip if we deliver from nophix and the parent dir is not nophix (or vice versa)
        pdir = os.path.basename(os.path.dirname(fqfile))
        if deliver_nophix and pdir != "nophix":
            continue
        if not deliver_nophix and pdir != run_name:
            continue
        
        # Skip if a compressed version of the current file exists
        if os.path.exists("{:s}.gz".format(fqfile)):
            print("WARNING: Both compressed and non-compressed versions of {:s} exists! " \
                  "Is compression/decompression in progress? Will deliver compressed version " \
                  "but you should make sure that the delivered files are complete!".format(fqfile))
            continue
        
        print("DEBUG: source_delivery_path = {:s}".format(os.path.dirname(fqfile)))
        
        fname = os.path.basename(fqfile)
        print(fname)
        # Split the file name according to CASAVA convention
        m = re.match(r'(\S+?)_[ACGTN\-]+_L0*(\d+)_R(\d)_\d+\.fastq(.*)', fname)
        if m is not None:
            lane = m.group(2)
            read = m.group(3)
            ext = m.group(4)
        else:
            # Split the file name according to bcbb convention
            m = re.match(r'(\d+)_(\d+)_([^_]+)_(\d+)_(?:nophix_)?(\d+)_fastq.txt(.*)', fname)
            if m is None:
                raise ValueError("Could not parse file name {:s} correctly! Exiting!".format(fname))
            lane = m.group(1)
            read = m.group(5)
            ext = m.group(6)
            
        dest_run_path = os.path.join(dest_proj_path, sample_name, run_name)
        dest_file_name = "{:s}.fastq{:s}".format("_".join([lane,
                                                           date,
                                                           fc_id,
                                                           sample_name,
                                                           read]),
                                                 ext)
        to_copy.append([fqfile,
                        dest_run_path,
                        dest_file_name])
    return to_copy

def rsync_files(to_copy, logfile, dry):
    # Iterate over the files to copy and create directories and copy files as necessary 
    successful = 0
    for src_file, dst_dir, dst_name in to_copy:
        dst_file = os.path.join(dst_dir, dst_name)
        print "Will copy (rsync) ", src_file, "to ", dst_file
        if not dry:
            
            # Create the destination directory if necessary
            logfile.write("[{:s}] - Creating run-level delivery directory: {:s} " \
                          "(or leaving it in place if already present)\n".format(utc_time(),
                                                                                 dst_dir))
            if os.path.exists(dst_dir):
                print("Directory {:s} already exists!".format(dst_dir))
            else:
                try:
                    # Create directory hierarchy with ug+rwX permissions
                    os.makedirs(dst_dir, 0770)
                except:
                    print("Could not create run-level delivery directory!")
                    sys.exit(1)
            
            # Rsync the file across 
            command_to_execute = ['rsync',
                                  '-ac',
                                  src_file,
                                  dst_file]
            
            logfile.write("[{:s}] - Executing command: {:s}\n".format(utc_time(), " ".join(command_to_execute)))
            logfile.flush()
            try:
                check_call(command_to_execute)
            except CalledProcessError, e:
                logfile.write("[{:s}] - rsync exited with exit code {:d}\n".format(utc_time(), e.returncode))
                raise e
            
            logfile.write("[{:s}] - rsync exited with exit code 0\n".format(utc_time()))
            successful += 1
        
            print("{:d} of {:d} files copied successfully".format(successful,len(to_copy)))
    
            # Modify the permissions to ug+rw
            os.chmod(dst_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP)
  

def main():
    parser = argparse.ArgumentParser(description="A script to help doing the deliveries, now using the Casava directory structure. " \
                                     "The user is asked to provide a project ID, a run name, and an UPPMAX project")

    parser.add_argument('-c', '--casava-path', action="store", dest="caspath", default='/proj/a2010002/nobackup/illumina/', 
                        help="Specify a path to a Casava directory manually")
    parser.add_argument('-l', '--log-path', action="store", dest="logpath", default='/proj/a2010002/private/delivery_logs', 
                        help="Specify a path to a log file")
    parser.add_argument('-i', '--interactive', action="store_true", dest="interactive", default=False, 
                        help="Interactively select samples to be delivered")
    parser.add_argument('-d', '--dry-run', action="store_true", dest="dry", default=False, 
                        help="Dry run: nothing will be done")
    parser.add_argument('-a', '--deliver-all-fcs', action="store_true", dest="deliver_all_fcs", default=False, 
                        help="rsync samples from all flow cells. Default is to only deliver from specified flowcell")
    parser.add_argument('-p', '--nophix', action="store_true", dest="deliver_nophix", default=False, 
                        help="Deliver fastq files from nophix subdirectory. Default is to deliver from run directory")
    parser.add_argument('project_name', action='store', help="Project name to deliver, e.g. J.Doe_10_01")
    parser.add_argument('flowcell_id', action='store', help="Flowcell id to deliver, e.g. 120824_BD1915ACXX")
    parser.add_argument('uppmax_id', action='store', help="UPPMAX project id to deliver to, e.g. b2012001")
    args = parser.parse_args()

    if not args.project_name in os.listdir(args.caspath): 
        print("Could not find project. Check directory listing:")
        for f in os.listdir(args.caspath): 
            print(f)
        sys.exit(0)

    fcid = args.flowcell_id
    fcid_comp = fcid.split('_')
    if len(fcid_comp) > 2:
        fcid = fcid_comp[0] + '_' + fcid_comp[-1]
        print("FCID format too long, trying {:s}".format(fcid))

    dt = datetime.now()
    time_str = "_".join([str(dt.year),
                         str(dt.month),
                         str(dt.day),
                         str(dt.hour),
                         str(dt.minute),
                         str(dt.second)])

    logfilename = os.path.join(os.path.normpath(args.logpath),"{:s}.log".format(time_str)) 
    if not args.dry:
        logfile = open(logfilename, "w")
    else:
        logfile = sys.stdout
         
    logfile.write("[{:s}] - Project to move files for:\n{:s}\n".format(utc_time(), args.project_name))
    logfile.flush()

    proj_base_dir = os.path.join(args.caspath, args.project_name)
    skip_list = []
    if args.interactive:
        for sample_dir in os.listdir(proj_base_dir):
            if not os.path.isdir(os.path.join(proj_base_dir,sample_dir)):
                continue
            if not query_yes_no("Deliver sample {:s}?".format(sample_dir), default="no"):
                skip_list.append(sample_dir)
    
    created_proj_dir_name = fixProjName(args.project_name)
    del_path_top = '/proj/' +  args.uppmax_id + "/INBOX/" + created_proj_dir_name 

    to_copy = get_file_copy_list(proj_base_dir,
                                 del_path_top,
                                 fcid,
                                 args.deliver_all_fcs,
                                 args.deliver_nophix,
                                 skip_list,
                                 logfile,
                                 args.dry)
    rsync_files(to_copy,
                logfile,
                args.dry)
        
    if not args.dry:
        logfile.close()
        


if __name__ == "__main__":
    main()
    