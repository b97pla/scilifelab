#!/usr/bin/env python
"""
Collect QC metrics and upload to statusdb

Usage:
  runQC_to_statusdb.py <flow cell dir> <statusdb> <run info YAML> [--dry_run --verbose]

"""

import os
import sys
import argparse
import couchdb

from bcbio.log import logger, setup_logging, version
from bcbio.log import logger2 as log
from bcbio.qc import FlowcellQCMetrics
from bcbio.qc.qcreport import _save_obj

def main(fc_dir,  statusdb, run_info_yaml):
    work_dir = os.getcwd()
    config = dict()
    config["log_dir"] = os.path.join(work_dir, "log")
    setup_logging(config)

    parts = fc_dir.split("_")
    fc_date, fc_name = parts[0], parts[3].rstrip("/")
    try:
        qc_obj = FlowcellQCMetrics(fc_date, fc_name, run_info_yaml, fc_dir, fc_dir)
    except:
        print "something wrong with %s" % fc_dir
        sys.exit()
    try:
        couch = couchdb.Server(url="http://%s" % statusdb)
        log.info("Connecting to server at %s succeeded" % statusdb)
    except:
        log.warn("Connecting to server at %s failed" % statusdb)
        sys.exit()
    db=couch['qc']
    for s in qc_obj.sample.keys():
        obj = qc_obj.sample[s]
        log.info("Saving sample %s" % obj.name())
        _save_obj(db, obj, statusdb)
    # Save flowcell object
    _save_obj(db, qc_obj, statusdb)

if __name__ == "__main__":
    usage = """
    runQC_to_statusdb.py post_process_config flowcellid runinfo
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("flowcelldir", type=str)
    parser.add_argument("statusdb", type=str)
    parser.add_argument("runinfo", type=str)

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                      default=False)
    parser.add_argument("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)

    args = parser.parse_args()

    main(args.flowcelldir, args.statusdb, args.runinfo)
    
