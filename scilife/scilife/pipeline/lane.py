"""Lane handling specific to SciLife"""

import os
import glob
from bcbio.solexa.flowcell import get_flowcell_info

def make_lane_items(dirs, config, run_items):
    """Generate lane_items list as in automated_initial_analysis"""
    lane_items = []
    (name, date) = get_flowcell_info(dirs['flowcell'])
    for info in run_items:
        dirname = os.path.join(dirs['flowcell'], "_".join([str(info['lane']), date, name]))
        sample_name = info.get("description", "")
        if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
            sample_name = "%s---%s" % (info.get("name", ""), sample_name)
        lane_name = "%s_%s_%s" % (info['lane'], date, name)
        if info.has_key('multiplex'):
            dirname = dirname + "_barcode"
            print "looking in " + dirname
            for bc in info['multiplex']:
                mname = bc['barcode_id']
                msample = bc['name']
                glob_str = "%s_%s_%s_%s*fastq.txt" % (info['lane'], date, name, mname)
                fq_files = glob.glob(os.path.join(dirname, glob_str))
                fastq1 = fq_files[0]
                fastq2 = fq_files[1] if fq_files[1] else None
                mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
                if msample is None:
                    msample = "%s---%s" % (sample_name, mname)
                lane_items.append((fastq1, fastq2, info['genome_build'], mname, msample, dirs, config ))
    return lane_items

def get_flowcell_id(run_info, fc_dir, check_bc=True, glob_ext="_fastq.txt"):
    lane = None
    for info in run_info:
        lane = info.get("lane", "")
    if check_bc:
        glob_str = "%s_*_barcode/*%s" % (lane, glob_ext)
    else:
        glob_str = "%s_*%s" % (lane, glob_ext)
    files = glob.glob(os.path.join(fc_dir, glob_str))
    try:
        (name, date) = get_flowcell_info(files[0])
    except:
        raise StandardError("No flowcell information found in " + str(fc_dir))
    return name, date
