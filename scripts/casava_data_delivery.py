# A script to help doing the deliveries.
# Now using the Casava directory structure
# The user is asked to provide a project ID, a run name, and an UPPMAX project

import sys
import os
import glob
import re
import grp
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
                 "_fastq.txt",
                 ".fastq..gz",
                 "_fastq.txt..gz"
                 ]
    for ext in fastq_ext:
        if fname.endswith(ext):
            return True
    return False

  
def create_final_name(fname, date, fc_id, sample_name):
    """Create the final name of the delivered file
    """
    
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
            raise ValueError("Could not parse file name {:s} correctly!".format(fname))
        lane = m.group(1)
        read = m.group(5)
        ext = m.group(6)
            
    dest_file_name = "{:s}.fastq{:s}".format("_".join([lane,
                                                       date,
                                                       fc_id,
                                                       sample_name,
                                                       read]),
                                             ext.replace('..','.'))
    return dest_file_name
      
def get_file_copy_list(proj_base_dir, dest_proj_path, fcid, deliver_all_fcs, deliver_nophix, skip_list):
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
      
        dest_run_path = os.path.join(dest_proj_path, sample_name, run_name)
        dest_file_name = create_final_name(fname,date,fc_id,sample_name)
        to_copy.append([fqfile,
                        dest_run_path,
                        dest_file_name])
    return to_copy

def rsync_files(to_copy, logfile, group, dry):
    # Iterate over the files to copy and create directories and copy files as necessary 
    successful = 0
    uid = os.getuid()
    gid = os.getgid()
    if group is not None and len(group) > 0:
        gid = grp.getgrnam(group).gr_gid
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
                    os.chown(dst_dir,uid,gid)
                except:
                    print("Could not create run-level delivery directory!")
                    clean_exit(1,logfile,dry)
            
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
            os.chown(dst_file,uid,gid)
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
    parser.add_argument('-g', '--group', action="store", dest="group", default="uppmax",
                        help="Group membership to set on copied files")
    parser.add_argument('project_name', action='store', help="Project name to deliver, e.g. J.Doe_10_01")
    parser.add_argument('flowcell_id', action='store', help="Flowcell id to deliver, e.g. 120824_BD1915ACXX")
    parser.add_argument('uppmax_id', action='store', help="UPPMAX project id to deliver to, e.g. b2012001")
    args = parser.parse_args()

    if not args.project_name in os.listdir(args.caspath): 
        print("Could not find project. Check directory listing:")
        for f in os.listdir(args.caspath): 
            print(f)
        clean_exit(0,None,args.dry)

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
                                 skip_list)
    
    # Prompt user if any of the files are non-compressed
    for fqfile, _, _ in to_copy:
        if os.path.splitext(fqfile)[1] == ".gz":
            continue
        print("WARNING: The file {:s}, which you are about to deliver, does not seem to be compressed. " \
              "It is recommended that you compress files prior to delivery.".format(fqfile))
        if query_yes_no("Do you wish to continue delivering " \
                        "uncompressed fastq files?", default="yes"):
            break
        clean_exit(1,logfile,args.dry)
            
    rsync_files(to_copy,
                logfile,
                args.group,
                args.dry)
        
    clean_exit(0,logfile,args.dry)
        
def clean_exit(exitcode, logfile, dry=False):
    """Close the logfile and exit with the given exit code
    """
    if not dry and logfile is not None:
        logfile.close()
    sys.exit(exitcode)

if __name__ == "__main__":
    main()
    
########## Tests ###########

import unittest
import shutil
import tempfile
import random
import uuid

class TestDataDelivery(unittest.TestCase):
    
    def test_fixProjName(self):
        """Fix project name
        """
        
        test_pnames = [("j.doe_11_01","J.Doe_11_01"),
                       ("j.Doe_11_01","J.Doe_11_01"),
                       ("J.doe_11_01","J.Doe_11_01"),
                       ("J.Doe_11_01","J.Doe_11_01"),
                       ("doe_11_01","Doe_11_01"),
                       ("j.d.doe_11_01","J.D.Doe_11_01"),]
        
        for test_pname, exp_pname in test_pnames:
            obs_pname = fixProjName(test_pname)
            self.assertEqual(obs_pname,
                             exp_pname,
                             "Did not get the expected fix ({:s}) for project name {:s} (got {:s})".format(exp_pname,test_pname,obs_pname))
        
    def test_is_fastq(self):
        """Determine if a file name corresponds to a fastq file
        """
        
        test_fnames = [("foo.fastq",True),
                       ("foo.fastq.gz",True),
                       ("foo_fastq.txt",True),
                       ("foo_fastq.txt.gz",True),
                       ("foo.fastq.bar",False),
                       ("foo.txt",False),]
        
        for test_fname, exp_result in test_fnames:
            obs_result = is_fastq(test_fname)
            self.assertEqual(obs_result,
                             exp_result,
                             "Did not get expected result ({:s}) for file name {:s}".format(str(exp_result),test_fname))
    
    def _create_test_files(self, root):
        
        to_copy = []
        for n in xrange(10):
            fd, sfile = tempfile.mkstemp(suffix=".tmp", prefix="rsync_test_", dir=root)
            os.close(fd)
            # Generate destination file hierarchies
            ddir = root
            for l in xrange(random.randint(1,5)):
                ddir = os.path.join(ddir,str(uuid.uuid4()))
            to_copy.append([sfile,ddir,"{:s}.tmp".format(str(uuid.uuid4()))])
        
        return to_copy
    
    def test_rsync_files(self):
        """Test the rsync functionality
        """
        root = tempfile.mkdtemp(prefix="rsync_test_")
        
        # Create some files to move
        to_copy = self._create_test_files(root)
        
        # Run rsync
        with open(os.devnull, 'w') as f:
            old_stdout = sys.stdout
            sys.stdout = f
            rsync_files(to_copy,sys.stdout,None,False)
            sys.stdout = old_stdout
            
        
        # Verify the copy process
        for src, ddir, dname in to_copy:
            self.assertTrue(os.path.exists(src),
                            "The rsync process have removed source file")
            self.assertTrue(os.path.exists(ddir) and os.path.isdir(ddir),
                            "The expected destination directory was not created")
            dfile = os.path.join(ddir,dname)
            self.assertTrue(os.path.exists(dfile) and os.path.isfile(dfile),
                            "The expected destination file was not created")
            exp_stat = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP
            obs_stat = stat.S_IMODE(os.stat(dfile).st_mode)
            self.assertEqual(obs_stat,
                             exp_stat,
                            "The mode of the created file is not as expected")
            
        shutil.rmtree(root)
        
    def test_rsync_set_group(self):
        """Test setting the group membership on rsync'd files
        """
        
        root = tempfile.mkdtemp(prefix="rsync_test_set_group_")
        avail_groups = os.getgroups()
        exp_group = grp.getgrgid(avail_groups[random.randint(1,len(avail_groups))-1])[0]
        
        # Create some files to move
        to_copy = self._create_test_files(root)
        
        # Run rsync
        with open(os.devnull, 'w') as f:
            old_stdout = sys.stdout
            sys.stdout = f
            rsync_files(to_copy,sys.stdout,exp_group,False)
            sys.stdout = old_stdout
            
        # Verify the copy process set the correct group on created directories
        for ddir in set([d[1] for d in to_copy]):
            gid = os.stat(ddir).st_gid
            obs_group = grp.getgrgid(gid)[0]
            self.assertEqual(obs_group,
                             exp_group,
                             "Failed to set group '{}' on directory. Group is {}".format(exp_group,
                                                                                    obs_group))
        
        # Verify the copy process set the correct group
        for src, ddir, dname in to_copy:
            dfile = os.path.join(ddir,dname)
            gid = os.stat(dfile).st_gid
            obs_group = grp.getgrgid(gid)[0]
            self.assertEqual(obs_group,
                             exp_group,
                             "Failed to set group '{}' on file. Group is {}".format(exp_group,
                                                                                    obs_group))
        
        
    def test_create_final_name(self):
        """Create the destination file name
        """
        
        date = "111111"
        fcid = "A11A22BCXX"
        sample_name = "P101_150B_index5"
        
        test_names = [("1_{}_{}_1_nophix_1_fastq.txt.gz".format(date,fcid),
                       "1_{}_{}_{}_1.fastq.gz".format(date,fcid,sample_name)),
                      ("1_{}_{}_1_nophix_1_fastq.txt".format(date,fcid),
                       "1_{}_{}_{}_1.fastq".format(date,fcid,sample_name)),
                      ("1_{}_{}_1_1_fastq.txt.gz".format(date,fcid),
                       "1_{}_{}_{}_1.fastq.gz".format(date,fcid,sample_name)),
                      ("{}_CGATGT_L001_R1_001.fastq.gz".format(sample_name),
                       "1_{}_{}_{}_1.fastq.gz".format(date,fcid,sample_name)),
                      ("{}_CGATGT_L001_R1_001.fastq..gz".format(sample_name),
                       "1_{}_{}_{}_1.fastq.gz".format(date,fcid,sample_name)),
                      ("{}_CGATGT_L001_R1_001.fastq".format(sample_name),
                       "1_{}_{}_{}_1.fastq".format(date,fcid,sample_name))]
        
        for test_fname, exp_result in test_names:
            obs_result = create_final_name(test_fname,date,fcid,sample_name)
            self.assertEqual(obs_result,
                             exp_result,
                             "Did not get expected final name ({:s}) for file name {:s}".format(exp_result,test_fname))
    
        # Try some illegal file names and assert that they raise exceptions
        test_names = ["1_{}_{}_1_nophix_1_fastq.gz".format(date,fcid),
                      "a_{}_{}_1_nophix_1_fastq.txt".format(date,fcid),
                      "{}_CGATRGT_L1_R1_001.fastq.gz".format(sample_name)]
        for test_name in test_names:
            with self.assertRaises(ValueError):
                create_final_name(test_name,date,fcid,sample_name)
            
    def test_get_file_copy_list(self):
        """Get list of files to copy and the destinations
        """
                     
        so = sys.stdout
        dn = open(os.devnull,"w")
        
        # Create a file hierarchy to search for files
        root = tempfile.mkdtemp(prefix="test_casava_data_delivery_")
        date = "111111"
        fcs = ["{}_{}".format(date,fcid) for fcid in ["FCA","FCB"]]
        
        # Create some sample files
        exp_files = []
        samples = []
        for n in xrange(2):
            sample = tempfile.mkdtemp(dir=root)
            samples.append(os.path.basename(sample))
            for fcid in fcs:
                fcdir = os.path.join(sample,fcid)
                nophixdir = os.path.join(fcdir,"nophix")
                for d in [fcdir,nophixdir]:
                    os.makedirs(d)
                    test_names = ["{:d}_{:s}_1_1_fastq.txt.gz".format(random.randint(1,8),
                                                                           fcid),
                                  "{}_CGATGT_L001_R1_001.fastq.gz".format(samples[-1]),
                                  "{}_CGATGT_L001_R1_001.fastq..gz".format(samples[-1]),]
                    for test_name in test_names:
                        test_file = os.path.join(d,test_name)
                        open(test_file,"w").close()
                        exp_files.append([samples[-1],
                                          fcid,
                                          os.path.basename(d) == "nophix",
                                          test_file,
                                          os.path.join(samples[-1],fcid),
                                          create_final_name(os.path.basename(test_name),date,fcid.split("_")[-1],samples[-1])])
                    
        # Get the list of files to copy under various conditions
        
        for deliver_all_fcs in [False, True]:
            for fcid in fcs:
                for deliver_nophix in [False, True]:
                    for skip_sample_list in [[],[samples[0]],[samples[1]],samples]:
                        sys.stdout = dn
                        obs_to_copy = sorted(get_file_copy_list(root,"",fcid,deliver_all_fcs,deliver_nophix,skip_sample_list))
                        sys.stdout = so
                        exp_to_copy = sorted([ef[3:6] for ef in exp_files if (deliver_all_fcs or ef[1] == fcid) and \
                                              deliver_nophix == ef[2] and \
                                              ef[0] not in skip_sample_list])
                        #import pdb; pdb.set_trace()
                        self.assertListEqual(obs_to_copy,
                                             exp_to_copy,
                                             "The files to copy result did not match the expected for " \
                                             "{:s}".format(", ".join(["{:s}: {:s}".format(k,v) for k, v in \
                                                                      dict(zip(["deliver_all_fcs",
                                                                                "fcid",
                                                                                "deliver_nophix",
                                                                                "skip_samples"],
                                                                               [str(deliver_all_fcs),
                                                                                fcid,
                                                                                str(deliver_nophix),
                                                                                " ".join(skip_sample_list)])).items()])))
          
                        
        
            