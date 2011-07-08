"""Lane handling specific to SciLife"""

import os
import glob
from bcbio.solexa.flowcell import get_flowcell_info
from scilife.pipeline import log
from scilife.pipeline.fastq import get_multiplex_items

def make_lane_items(info, fc_date, fc_name, dirs, config):
    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", "")
    log.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    lane_items = []
    if multiplex:
        log.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
        mitems = get_multiplex_items(multiplex, info['lane'], dirs['fc_dir'], fc_name, fc_date)
        for fastq1, fastq2, mlane_name, msample in mitems:
            lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample, dirs, config))
    else:
        # TODO: Not multiplex: what to do?
        pass
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
