# A script to help doing the deliveries.
# Now using the Casava directory structure
# The user is asked to provide a project ID, a run name, and an UPPMAX project

import sys, os, yaml, glob, shutil
from datetime import datetime
import optparse

def_casava_path = '/proj/a2010002/nobackup/illumina/'
def_log_path = '/bubo/home/h9/mikaelh/delivery_logs/'

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

if len(sys.argv) < 4:
    print "USAGE: python " + sys.argv[0] + " <project ID> <flow cell ID, e g 120824_BD1915ACXX> <UPPMAX project> [-d Dry run -i Interactive -c <path to Casava dir> (optional)] [-l <path to log file dir [optional]>]"
    sys.exit(0)

parser = optparse.OptionParser()
parser.add_option('-c', '--casava-path', action="store", dest="caspath", default=def_casava_path, help="Specify a path to a Casava directory manually")
parser.add_option('-l', '--log-path', action="store", dest="logpath", default=def_log_path, help="Specify a path to a log file")
parser.add_option('-i', '--interactive', action="store_true", dest="interactive", default=False, help="Interactively select samples to be delivered")
parser.add_option('-d', '--dry-run', action="store_true", dest="dry", default=False, help="Dry run")
parser.add_option('-a', '--deliver-all-fcs', action="store_true", dest="deliver_all_fcs", default=False, help="rsync samples from all flow cells")

(opts, args) = parser.parse_args()

base_path = opts.caspath
log_path = opts.logpath
interactive = opts.interactive
dry = opts.dry
deliver_all_fcs = opts.deliver_all_fcs

# print "DEBUG: Base path: ", base_path

#projid = sys.argv[1].lower()
projid = sys.argv[1]
fcid = sys.argv[2]
uppmaxproj = sys.argv[3]

dt = datetime.now()
time_str  = str(dt.year) + "_" + str(dt.month) + "_" + str(dt.day) + "_" + str(dt.hour) + "_" + str(dt.minute) + "_" + str(dt.second)

if not dry:
    logfilename = log_path + time_str + ".log" 
    logfile = open(logfilename, "w")

print "Project to move files for:", projid
#if not dry: logfile.write("Project to move files for:" + "\n" + fixProjName(projid) + "\n")

if not projid in os.listdir(base_path): 
    print "Could not find project. Check directory listing:"
    for f in os.listdir(base_path): print f
    sys.exit(0)

if not dry: logfile.write("Project to move files for:" + "\n" + projid + "\n")

dirs_to_copy_from = []

#proj_base_dir = os.path.join(base_path, fixProjName(projid))
proj_base_dir = os.path.join(base_path, projid)

os.chdir(proj_base_dir)
run_found = False
for i in glob.glob('*'):
    if os.path.isdir(i):
        dirs_to_copy_from.append(i)
        
if not dry: logfile.flush()

# Create directory in user's INBOX

created_proj_dir_name = fixProjName(projid)
del_path_top = '/proj/' +  uppmaxproj + "/INBOX/" + created_proj_dir_name 

print "Will create a top-level project directory", del_path_top       
if not dry: 
    logfile.write("Creating top-level delivery directory:" + del_path_top + " (or leaving it in place if already present)\n")
    if os.path.exists(del_path_top):
        print "Directory", del_path_top, " already exists!"
    else:
        try:
            os.mkdir(del_path_top)
        except:
            print "Could not create project-level delivery directory!"
            sys.exit(0)

skip_list = []
if interactive:
    for sample_dir in dirs_to_copy_from:
        print "Deliver sample ", sample_dir, " ? (y/n)"
        ans = raw_input()
        if ans.lower() != 'y': skip_list.append(sample_dir)
   
for sample_dir in dirs_to_copy_from:
    if sample_dir in skip_list: continue
    sample_path = os.path.join(del_path_top, sample_dir)
    if not dry: 
        logfile.write("Creating sample-level delivery directory:" + sample_path + " (or leaving it in place if already present)\n")
        if os.path.exists(sample_path):
            print "Directory ", sample_path, " already exists!"
        else:
            try:
                os.mkdir(sample_path)
            except:
                print "Could not create sample-level delivery directory!"
                sys.exit(0)

    if not deliver_all_fcs:            
        phixfiltered_path = os.path.join(proj_base_dir, sample_dir, fcid, "nophix")
    else:
        phixfiltered_path = os.path.join(proj_base_dir, sample_dir, "*", "nophix")
    for fq in glob.glob(os.path.join(phixfiltered_path, "*fastq*")):
        [path, fname] = os.path.split(fq)
        run_dir = os.path.split(os.path.split(fq)[0])[0]
        # print "DEBUG: Run dir = ", run_dir
        run_name = os.path.basename(os.path.split(os.path.split(fq)[0])[0])
        if not dry:
            if not os.path.exists(os.path.join(sample_path, run_name)):
                try:
                    os.mkdir(os.path.join(sample_path, run_name))
                except:
                    print "Could not create run level directory!"
                    print os.path.join(sample_path, run_name)
                    sys.exit(0)
        sample = os.path.basename(sample_path)
        print fname
        [lane, date, fcid, bcbb_id, nophix, read, dummy] = fname.split('_') # e.g. 4_120821_BC118PACXX_1_nophix_2_fastq.txt
        dest_file_name = lane + "_" + date + "_" + fcid + "_" + sample + "_" + read + ".fastq" 
        dest = os.path.join(sample_path, run_name, dest_file_name)
        print "Will copy (rsync) ", fq, "to ", dest 
        if not dry: 
            command_to_execute = 'rsync -ac ' + fq + ' ' + dest
            os.system(command_to_execute) 
            logfile.write("Executing command: " + command_to_execute)
            logfile.flush()

if not dry: 
    os.chdir(del_path_top)
    logfile.close()
    os.system("chmod -R g+rw " + del_path_top)
