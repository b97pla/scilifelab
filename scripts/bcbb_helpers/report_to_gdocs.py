#!/usr/bin/env python
"""Helper script that uploads demultiplex (barcode) counts and QC data to godcs.
"""
import os
import glob
import yaml
import sys

from optparse import OptionParser
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.demultiplex import _find_demultiplex_stats_htm
from bcbio.pipeline.demultiplex import _parse_demultiplex_stats_htm
from bcbio.pipeline.demultiplex import _write_demultiplex_metrics
from bcbio.google.sequencing_report import create_report_on_gdocs


def main(run_id, config_file, run_info_file=None, dryrun=False):
    assert os.path.exists(config_file), \
    "The post process configuration file, %s, could not be found" % config_file

    config = load_config(config_file)
    archive_dir = config["analysis"]["store_dir"]
    archive_dir = os.path.normpath(archive_dir)
    analysis_dir = config["analysis"]["base_dir"]
    analysis_dir = os.path.normpath(analysis_dir)

    dirs = {"work": os.path.join(analysis_dir, run_id),
            "flowcell": os.path.join(archive_dir, run_id)}

    if run_info_file is None:
        run_info_file = os.path.join(dirs["flowcell"], "run_info.yaml")

    assert run_id, \
    "No run id was specified"
    assert "gdocs_upload" in config, \
    "The configuration file, %s, has no section specifying the Google docs details" % config_file
    assert os.path.exists(run_info_file), \
    "The run info configuration file, %s, could not be found" % run_info_file
    assert os.path.exists(dirs["flowcell"]), \
    "The flowcell directory, %s, could not be found" % dirs["flowcell"]
    assert os.path.exists(dirs["work"]), \
    "The work directory, %s, could not be found" % dirs["work"]

    fc_name, fc_date = get_flowcell_info(dirs["flowcell"])
    # If we have no bc_metrics files in the workdir, we may be looking at a Casava run.
    # In that case, attempt to parse the Demultiplex_Stats.htm file and create bc_metrics files
    bc_metrics_files = glob.glob(os.path.join(dirs["work"], "*_barcode", "*bc.metrics"))
    bc_metrics_files += glob.glob(os.path.join(dirs["work"], "*bc.metrics"))
    if len(bc_metrics_files) == 0:
        casava_stats = _find_demultiplex_stats_htm(run_id, config)
    else:
        casava_stats = None

    info = []
    if casava_stats:
        bc_metrics = _parse_demultiplex_stats_htm(casava_stats)
        with open(run_info_file) as fh:
            info = yaml.load(fh)

    for item in info:
        metrics_file = "_".join([item["lane"], fc_date, fc_name, "bc.metrics"])
        multiplex = item.get("multiplex", [])
        for plex in multiplex:
            plex["lane"] = item["lane"]

        _write_demultiplex_metrics(multiplex, bc_metrics, os.path.join(dirs["work"], metrics_file))

    print("A report will be created on Google Docs based on the demultiplexed data in {}".format(dirs["work"]))
    print("The configuration file is {0} and the run info file is {1}".format(config_file, run_info_file))
    print("The run was started on {0} and has flowcell id {1}".format(fc_date, fc_name))

    if not dryrun:
        create_report_on_gdocs(fc_date, fc_name, run_info_file, dirs, config)
    else:
        print("DRY-RUN: nothing uploaded")


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-c", "--config-file", dest="config_file", default=None)
    parser.add_option("-f", "--run-info-file", dest="run_info_file", default=None)
    parser.add_option("-n", "--dry-run", dest="dryrun", action="store_true", default=False)
    options, args = parser.parse_args()

    run_id = None
    if len(args) == 1:
        run_id = args[0]
    else:
        print(__doc__)
        sys.exit()

    main(run_id, options.config_file,
         os.path.normpath(options.archive_dir), os.path.normpath(options.analysis_dir),
         options.run_info_file, options.dryrun)
