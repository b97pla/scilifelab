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

def map_fastq_barcode(dirs, run_items):
    """Data has been delivered to fastq dir, in which barcode ids have
    been substituted with barcode names as defined in YAML run information.
    Here files in fastq_dir are mapped to their original barcode names"""
    glob_str = "*.fastq"
    files = glob.glob(os.path.join(dirs['fastq'], glob_str))
    fc_dir = {}
    (name, date) = get_flowcell_info(files[0])
    fc = "_".join([date, name])
    if not fc_dir.has_key(fc):
        fc_dir[fc] = fc
    if len(fc_dir.keys()) > 1:
        raise StandardError("Can currently only process one flowcell per directory")
    dirs['flowcell'] = os.path.join(dirs['align'], os.pardir, fc_dir.keys()[0])
    if not os.path.exists(dirs['flowcell']):
        os.makedirs(dirs['flowcell'])
    d = dict()
    for item in run_items:
        d[item['lane']] = item
    lane = None
    lanes = dict()
    barcodeids = dict()
    sample = None
    src2tgt = dict()
    mplexinfo = collections.defaultdict(list)
    for f in files:
        src = os.path.basename(f)
        parts = src.split("_")
        lane = int(parts[0])
        

        lanes[lane] = lane
        if not mplexinfo.has_key(lane):
            mplexinfo[lane] = list()
            barcodeids[lane] = dict()
        lane_info = None
        if d[lane].has_key('multiplex'):
            dirname = os.path.join(dirs['flowcell'], "%s_%s_%s_barcode" % (lane, date, name))
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            mplex = d[lane]['multiplex']
            for mp in mplex:
                if src.find(mp['name']) != -1:
                    tgt = src.replace(mp['name'], mp['barcode_id'])
                    tgt = tgt.replace('.fastq', '_fastq.txt')
                    os.symlink(os.path.join(dirs['fastq'], src), os.path.join(dirname, tgt))
                    if not barcodeids[lane].has_key(mp['barcode_id']):
                        mplexinfo[lane].append(mp)
                    barcodeids[lane][mp['barcode_id']] = True
        else:
            pass
    new_run_items = list()
    for l in lanes.keys():
        if d.has_key(l):
            d[l]['multiplex'] = mplexinfo[l]
            new_run_items.append(d[l])
    return (dirs, new_run_items)

def organize_fastq_samples():
    """Organize delivered fastq files by sample name, handling multiplexing"""
    pass


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


