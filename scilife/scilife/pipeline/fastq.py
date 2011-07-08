"""Pipeline utilities

- get_fastq_files
  Modified from bcbio since our regexp is *.fastq

- get_samples_from_fastq_dir
  Based on fastq names formatted as lane_date_flowcell_sample_[indexI]_1/2.fastq,
  retrieves the name, date, and sample
  USE: if run_info not present?

"""
import os
import glob 
import subprocess
import collections

from bcbio.solexa.flowcell import get_flowcell_info

def get_fastq_files(directory, sample, fc_name, bc_name=None):
    """Retrieve fastq files for the given sample, ready to process.
    """
    if bc_name:
        glob_str = "*_*_%s_%s_%s_*.fastq*" % (fc_name, sample, bc_name)
    else:
        glob_str = "*_*_%s_%s_*.fastq*" % (fc_name, sample)
    files = glob.glob(os.path.join(directory, glob_str))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, sample, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

def get_samples_from_fastq_dir(fastq_dir, info=None):
    """Get information about samples in a directory with fastq files.

    The function will search for files that comply with given 
    regular expressions. Therefore, it will only work given that certain
    naming conventions are fulfilled. For instance, it assumes that 
    all fastq files are from the same flowcell.

    The function returns a flowcell flowcell name, flowcell date, and 
    """
    run_items = list()
    samples = {}
    if info:
        pass
    else:
        fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq*"))
        if not fastq_files:
            raise StandardError("No fastq files found in " + fastq_dir)
        fcnames = {}
        fcdates = {}
        for fq in fastq_files:
            (name, date) = get_flowcell_info(fastq_files[0])
            fcnames[name] = name
            fcdates[date] = date
        if len(fcnames.keys()) > 1:
            raise StandardError("A directory may only contain files from one flowcell")
    for fq in fastq_files:
        (_, fc) = os.path.split(fq)
        # files have to be lane_date_flowcell_sample_1/2.fastq
        sample = fc.split("_")[3]
        if not samples.has_key(sample):
            run_items.append({'sample':sample})
        samples[sample] = True
    return name, date, run_items

# TODO: these two could probably also be handled much more efficiently
def get_barcoded_project_files(multiplex, lane, fc_dir, fc_name):
    fq = list()
    for bc in multiplex:
        fq.append(_get_fastq_files(fc_dir, lane, fc_name, ".fastq", bc['name']))
    return fq

def get_barcoded_fastq_files(multiplex, lane, fc_dir, fc_name, fc_date):
    fq = list()
    bc_dir = "%s_%s_%s_barcode" % (lane, fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    for bc in multiplex:
        if not os.path.exists(bc_dir):
            raise IOError("No barcode directory found: " + str(bc_dir))
        fq.append(_get_fastq_files(bc_dir, lane, fc_name, bc_name=bc['barcode_id']))
    return fq

def get_multiplex_items(multiplex, lane, fc_dir, fc_name, fc_date):
    mitems = list()
    bc_dir = "%s_%s_%s_barcode" % (lane, fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    lane_name = "%s_%s_%s" % (lane, fc_date, fc_name)
    for bc in multiplex:
        mname = bc['barcode_id']
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        msample = bc['name']
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        if not os.path.exists(bc_dir):
            raise IOError("No barcode directory found: " + str(bc_dir))
        fastq1, fastq2 = _get_fastq_files(bc_dir, lane, fc_name, bc_name=bc['barcode_id'])
        mitems.append((fastq1, fastq2 , mlane_name, msample))
    return mitems
    

def get_single_fastq_files(lane, fc_dir, fc_name):
    return _get_fastq_files(fc_dir, lane, fc_name)

# Note: do I even need to use this? The only difference from bcbio.pipeline.fastq is the first glob_str
# There, it seems as if one always has paired end reads?!?
def _get_fastq_files(directory, lane, fc_name, glob_ext="_fastq.txt", bc_name=None):
    """Retrieve processed fastq files, possibly demultiplexed"""
    if bc_name:
        glob_str = "%s_*%s_%s_*%s" % (lane, fc_name, bc_name, glob_ext)
    else:
        glob_str = "%s_*%s*_%s" % (lane, fc_name, glob_ext)
    files = glob.glob(os.path.join(directory, glob_str))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)



# TODO: these two could probably be handled much more efficiently
def convert_barcode_id_to_name(multiplex, fc_name, fq):
    """Convert barcode id to sample description, changing extension from _fastq.txt to .fastq in the process"""
    fqout = list([None, None])
    bcid2name = dict([(mp['barcode_id'], mp['name']) for mp in multiplex])
    for bcid in bcid2name.keys():
        mstr = "%s_%s_" % (fc_name, bcid) 
        if fq[0].find(mstr) != -1:
            from_str = "%s_%s_" %(fc_name, bcid)
            to_str   = "%s_%s_" %(fc_name, bcid2name[bcid])
            fqout[0] = fq[0].replace(from_str, to_str)
            if not fq[1] == None:
                fqout[1] = fq[1].replace(from_str, to_str)
    fqout[0] = fqout[0].replace("_fastq.txt", ".fastq")
    if not fqout[1] == None:
        fqout[1] = fqout[1].replace("_fastq.txt", ".fastq")
    return os.path.basename(fqout[0]), (os.path.basename(fqout[1]) if len(fqout) > 1 else None)

def convert_name_to_barcode_id(multiplex, fc_name, fq):
    """Convert sample description to barcode id, changing extension from .fastq to _fastq.txt in the process"""
    fqout = list([None, None])
    name2bcid = dict([(mp['name'], mp['barcode_id']) for mp in multiplex])
    for name in name2bcid.keys():
        mstr = "%s_%s_" % (fc_name, name) 
        if fq[0].find(mstr) != -1:
            from_str = "%s_%s_" %(fc_name, name)
            to_str   = "%s_%s_" %(fc_name, name2bcid[name])
            fqout[0] = fq[0].replace(from_str, to_str)
            if not fq[1] == None:
                fqout[1] = fq[1].replace(from_str, to_str)
    fqout[0] = fqout[0].replace(".fastq", "_fastq.txt")
    if not fqout[1] == None:
        fqout[1] = fqout[1].replace(".fastq", "_fastq.txt")
    return os.path.basename(fqout[0]), (os.path.basename(fqout[1]) if len(fqout) > 1 else None)

