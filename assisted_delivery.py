# A script to help doing the deliveries.

# The user is asked to provide a project ID, a run name, and an UPPMAX project

import sys, os, yaml, glob, shutil
from datetime import datetime

if len(sys.argv) < 5:
    print "USAGE: python " + sys.argv[0] + " <project ID> <run name> <UPPMAX project> <Interactive, y/n>"
    sys.exit(0)

interactive = True

projid = sys.argv[1].lower()
runname = sys.argv[2].strip("/")
yamlfile = 'store/' + runname + "/run_info.yaml"
uppmaxproj = sys.argv[3]

if sys.argv[4].lower() == "n": interactive = False

projdata = yaml.load(open(yamlfile))

dt = datetime.now()
time_str  = str(dt.year) + "_" + str(dt.month) + "_" + str(dt.day) + "_" + str(dt.hour) + "_" + str(dt.minute) + "_" + str(dt.second)
logfilename = "/bubo/home/h9/mikaelh/delivery_logs/" + time_str + ".log" 
logfile = open(logfilename, "w")

print "Project to copy files for:", projid
logfile.write("Project to copy files for:" + "\n" + projid + "\n")

matching = []
skipped = []

for entry in projdata:

    if entry['description'].split(',')[-1].strip().lower()==projid: 
        matching.append(entry['lane'])
    else:
        skipped.append(entry['lane'])

if interactive:
    print "I will copy files from lanes " + ",".join(matching) + " and skip lanes " + ",".join(skipped)
    logfile.write("Files will be copied from lanes" + ",".join(matching) + "\n")
    if interactive:
        r = raw_input("Are you happy with this? (y/n)")
        if r.lower() == "n": 
            print "Exiting ..."
            logfile.write("Aborted by user")
            logfile.close()
            sys.exit(0)

logfile.flush()
# Create directory in user's INBOX

temp = runname.split('_')
start_date = temp[0]
flow_cell = temp[3][0] # A or B
created_dir_name =  "20" + start_date + flow_cell + "_hiseq2000"

del_path = '/bubo/proj/' +  uppmaxproj + "/INBOX/" + created_dir_name

print "I will now create a delivery directory", del_path       
logfile.write("Creating delivery directory:" + del_path + " (or leaving it in place if already present)\n")

if interactive:
    r = raw_input("Type n to abort... ")
    if r.lower() == "n": 
        print "Aborted"
        logfile.write("Aborted by user")
        logfile.close()
        sys.exit(0)

if os.path.exists(del_path):
        print "Directory", del_path, " already exists!"
else:
    try:
        os.mkdir(del_path)
    except:
        print "Could not create delivery directory!"
        sys.exit(0)

# Start looking for the files to transfer

temp = runname.split('_')
dirs_to_process = []
for m in matching:
    d = m + "_" + temp[0] + "_" + temp[3] 
    dirs_to_process.append(d)

os.chdir(runname)

for d in dirs_to_process:
    dirpath = d + "_barcode"
    print "Examining directory", os.getcwd() + "/" + dirpath
    if not os.path.exists(dirpath):
        print "Could not find directory", dirpath 
        sys.exit(0)
    os.chdir(dirpath)
    bcname = d + "_bc.metrics"
    lane = dirpath[0]
    print "LANE ", lane
    logfile.write("LANE " + lane + "\n")
# Print table of Illumina vs. bcbb barcodes
    sample_id_and_idx = {}
    lane_info = "none"
    for entry in projdata:
        if entry['lane'] == lane:
            lane_info = entry
    
    if lane_info.has_key('multiplex'):
        for bc in lane_info['multiplex']:
            sample_id_and_idx[bc['barcode_id']] = bc['name']
        print sample_id_and_idx

        print "Pipeline index\tIllumina index/sample ID\tMatches"
        logfile.write("Pipeline index\tIllumina index/sample ID\tMatches\n")
        if os.path.exists(bcname):
            for line in open(bcname):
                [bcbb_bc, hits] = line.strip().split()
                try:
                    print bcbb_bc + "\t" + sample_id_and_idx[int(bcbb_bc)] + "\t" + hits
                    logfile.write(bcbb_bc + "\t" + sample_id_and_idx[int(bcbb_bc)] + "\t" + hits + "\n")
                except:
                    if bcbb_bc == "unmatched": 
                        print bcbb_bc + "\t" + "N.A." + "\t" + hits 
                        logfile.write(bcbb_bc + "\t" + "N.A." + "\t" + hits + "\n")
                    else:
                        print "Encountered parsing error in barcode conversion: " + bcbb_bc
                        print sample_id_and_idx
                        sys.exit(0)              
        else:
            print "BC metrics file", bcname, " not found"
            sys.exit(0)

    else:
        print "Non-multiplexed lane"
        logfile.write("Non-multiplexed lane\n")

    # print os.listdir(".")
    files_to_copy = []

    for fastq_file in glob.glob("*fastq.txt"):
        if lane_info.has_key('multiplex'):
            if 'unmatched' in fastq_file: continue
            # Extract barcode
            [lane, date, run_id, bcbb_bc, pe_read, dummy] = fastq_file.split("_")
            customer_sample_id = sample_id_and_idx[int(bcbb_bc)]
            new_file_name = lane + "_" + date + "_" + run_id + "_" + customer_sample_id.replace("/", "_") + "_" + pe_read + ".fastq"   
        else:
            # 8_110707_AB02J6ABXX_trim_1_fastq.txt
            
            [lane, date, run_id, name, pe_read,dummy] = fastq_file.split("_")
            new_file_name = lane + "_" + date + "_" + run_id + "_" + name + "_" + pe_read + ".fastq"   
        print "Preparing to copy file", fastq_file, "as ", new_file_name
        files_to_copy.append([fastq_file, new_file_name])

    if interactive:
        r = raw_input("Go ahead with copying? (press n to abort)")
        if r.lower() == "n": 
            logfile.write("Aborted by user")
            logfile.close()
            sys.exit(0)
    
    for pair in files_to_copy:
        source = os.getcwd() + "/" + pair[0]
        dest = del_path + "/" + pair[1]
        print "Copying ", source, "to", dest
        logfile.write("Copying " + source + " to " + dest + "\n")
        logfile.flush()
        shutil.copyfile(source, dest)

    os.chdir("..")
    if interactive: raw_input("Press a key to go to the next lane ...")

logfile.close()
