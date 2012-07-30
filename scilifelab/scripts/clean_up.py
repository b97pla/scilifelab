#!/usr/bin/env python

import sys, os, yaml, glob, shutil
from datetime import datetime

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


def main(projid, runname, abbr_runname, uppmaxproj, dry):

    base_path = '/proj/a2010002/nobackup/illumina/'
    base_yaml_path = '/proj/a2010002/archive/'

    yamlfile = base_yaml_path + runname + "/run_info.yaml"
    projdata = yaml.load(open(yamlfile))

    matching = set()
    available = set()


    # Builds a set containing the projects that match the run_info.yaml
    # description field against the supplied flowcell directory (run name)

    for entry in projdata:
        proj_name = entry['description'].split(',')[-1].strip()
        available.add(proj_name)
        if proj_name.lower()==projid:
            matching.add(entry['lane'])
        elif entry.has_key('multiplex'):
                for sample in entry['multiplex']:
                    if sample.has_key('sample_prj'):
                        available.add(sample['sample_prj'])
                        if sample['sample_prj'].split(',')[-1].strip().lower()==projid:
                            matching.add(entry['lane'])

    if len(matching)==0:
        print "No matching project found. Possibilities:"
        for prid in sorted(available):
            print prid
        sys.exit(0)

    temp = runname.split('_')
    dirs_to_process = []
    for match in sorted(matching):
        dir_nophix = match + "_" + temp[0] + "_" + temp[3] + "_nophix"
        dirs_to_process.append(dir_nophix)

    ## XXX
    os.chdir(base_path + runname)


    for d in dirs_to_process:
        dirpath = d + "_barcode"
        if not os.path.exists(dirpath):
            print "Could not find directory", dirpath
            print "Standing in ", os.getcwd()
            sys.exit(0)

        # Enter the demultiplexed lane directory (_barcode) and read barcode_id's
        os.chdir(dirpath)
        bcname = d + "_bc.metrics"
        lane = dirpath[0]
        lane_info = "none"
        for entry in projdata:
            if entry['lane'] == lane:
                lane_info = entry

        # Maps bcbio-nextgen internal barcode identifiers with sample names
        sample_id_and_idx = {}
        if lane_info.has_key('multiplex'):
            for bc in lane_info['multiplex']:
                if bc.has_key('sample_prj'):
                    if bc['sample_prj'].split(',')[-1].strip().lower() == projid:
                        sample_id_and_idx[bc['barcode_id']] = bc['name']

            # Check that demultiplexing has been done
            if not os.path.exists(bcname):
                print "BC metrics file", bcname, " not found"
                sys.exit(0)


        files_to_copy = []
        for fastq_file in glob.glob("*fastq.txt"):
            if lane_info.has_key('multiplex'):
                new_file_name = None
                if 'unmatched' in fastq_file:
                    continue

                [lane, date, run_id, nophix, bcbb_bc, pe_read, dummy] = fastq_file.split("_")

                if int(bcbb_bc) in sample_id_and_idx:
                    customer_sample_id = sample_id_and_idx[int(bcbb_bc)]
                    new_file_name = lane + "_" + date + "_" + run_id + "_" + customer_sample_id.replace("/", "_") + "_" + pe_read + ".fastq"
            else:
                [lane, date, run_id, nophix, name, pe_read,dummy] = fastq_file.split("_")
                new_file_name = lane + "_" + date + "_" + run_id + "__" + pe_read + ".fastq"

            if new_file_name != None:
                files_to_copy.append([fastq_file, new_file_name])

        for fastq_file in glob.glob("*fastq.txt"):
            key = fastq_file.split('_1_fastq.txt')[0].split('_2_fastq.txt')[0]
            to_delete = [os.getcwd() + "/" + fastq_file]
            to_delete.extend(glob.glob(base_path +runname+"/"+key+"*"))
            to_delete.extend(glob.glob(base_path +runname+"/fastqc/"+key+"*"))
            to_delete.extend(glob.glob(base_path +runname+"/fastq_screen/"+key+"*"))

            if dry:
                print "WILL DELETE"
                for f in to_delete:
                    if os.path.isfile(f):
                        print f
            else:
                print "DELETING"
                for f in to_delete:
                    if os.path.isfile(f):
                        os.remove(f)
                    elif os.path.isdir(f):
                        print "skipping... %s" % f

        os.chdir(base_path + runname)


if __name__ == "__main__":

    if len(sys.argv) < 5:
        print "USAGE: python " + sys.argv[0] + " <project ID> <run name> <UPPMAX project> <Dry run, y/n>"
        #XXX
        # Example print "USAGE: python " + sys.argv[0] + " t.olsson <run name> <UPPMAX project> <Dry run, y/n>"
        sys.exit(0)
    
    dry = True

    projid = sys.argv[1].lower()
    runname = sys.argv[2].strip("/")
    abbr_runname = runname.split("_")[0] + "_" + runname.split("_")[3]
    uppmaxproj = sys.argv[3]

    if sys.argv[4].lower() == "n":
        dry = False
    
    main(projid, runname, abbr_runname, uppmaxproj, dry)
