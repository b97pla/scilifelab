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

if len(sys.argv) < 5:
    print "USAGE: python " + sys.argv[0] + " <project ID> <run name> <UPPMAX project> <Dry run, y/n>  [-c <path to Casava dir> (optional)] [-l <path to log file dir [optional]>]"
    sys.exit(0)

parser = optparse.OptionParser()
parser.add_option('-c', '--casava-path', action="store", dest="caspath", default=def_casava_path, help="Specify a path to a Casava directory manually")
parser.add_option('-l', '--log-path', action="store", dest="logpath", default=def_log_path, help="Specify a path to a log file")

(opts, args) = parser.parse_args()

base_path = opts.caspath
log_path = opts.logpath

print "DEBUG: Base path: ", base_path

dry = True

#projid = sys.argv[1].lower()
projid = sys.argv[1]
runname = sys.argv[2].strip("/")
runname_comps = runname.split("_")
if len(runname_comps) == 4: abbr_runname = runname.split("_")[0] + "_" + runname.split("_")[3] 
else: abbr_runname = runname
uppmaxproj = sys.argv[3]

if sys.argv[4].lower() == "n": dry = False

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

if not dry: logfile.write("Project to move files for:" + "\n" + projid + "\n")

dirs_to_copy_from = []

#proj_base_dir = os.path.join(base_path, fixProjName(projid))
proj_base_dir = os.path.join(base_path, projid)

os.chdir(proj_base_dir)
for i in glob.glob('*'):
    run_found = False
    avail = set()
    # print "DEBUG: Abbreviated run name: ", abbr_runname
    if os.path.isdir(i):
        dir_cont = os.listdir(i)
        for d in dir_cont: avail.add(d) 
        if abbr_runname in dir_cont: 
            print "DEBUG: Sample ", i, " has been run on FC ", runname
            run_found = True
            dirs_to_copy_from.append(i)
        else:
            print "DEBUG: Sample ", i, " has not been run on FC ", runname

if not run_found:
    print "Could not find specified run for the samples in this project. Available runs:"
    for a in avail: print a
    sys.exit(0)

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

print "Will create sample-level directories"

for direc in dirs_to_copy_from:
    sample_path = os.path.join(del_path_top, direc)
    run_path = os.path.join(sample_path, abbr_runname)
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
        logfile.write("Creating run-level delivery directory:" + run_path + " (or leaving it in place if already present)\n")
        if os.path.exists(run_path):
            print "Directory ", run_path, " already exists!"  
        else:    
            try:
                os.mkdir(run_path)
            except:
                sys.exit("Could not create run-level delivery directory!")

    run_path = os.path.join(sample_path, abbr_runname)
    phixfiltered_path = os.path.join(proj_base_dir, direc, abbr_runname, "nophix")
    for fq in glob.glob(os.path.join(proj_base_dir, direc, abbr_runname, "nophix", "*fastq.txt")):
        [path, fname] = os.path.split(fq)
        dest_file_name = fname.replace("_fastq.txt", ".fastq")
        dest = os.path.join(run_path, dest_file_name)
        print "Will copy (rsync) ", fq, "to ", dest 
        if not dry: 
            command_to_execute = 'rsync -a ' + fq + ' ' + dest
            os.system(command_to_execute) 
            logfile.write("Executing command: " + command_to_execute)
            logfile.flush()

if not dry: 
    os.chdir(del_path_top)
    logfile.close()
    os.system("chmod -R g+rw " + del_path_top)
