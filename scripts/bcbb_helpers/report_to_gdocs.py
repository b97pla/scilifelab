#!/usr/bin/env python
"""Helper script that uploads demultiplex (barcode) counts and QC data to godcs.
"""
import os
import glob
import yaml
import sys
import tempfile
from collections import defaultdict

from optparse import OptionParser
from bcbio.pipeline.config_loader import load_config
import bcbio.pipeline.demultiplex as dmx 
import bcbio.solexa.flowcell as fc 
import bcbio.solexa.samplesheet as ssheet
from bcbio.google.sequencing_report import create_report_on_gdocs


def main(run_id, config_file, run_info_file=None, dryrun=False):
    
    assert run_id, \
    "No run id was specified"
    assert os.path.exists(config_file), \
    "The configuration file, {}, could not be found".format(config_file)

    config = load_config(config_file)    
    assert "gdocs_upload" in config, \
    "The configuration file, {}, has no section specifying the Google docs details".format(config_file)

    analysis_cfg = config.get("analysis",{})
    if "store_dir" in analysis_cfg:    
        archive_dir = os.path.join(analysis_cfg["store_dir"], run_id)
    else:
        archive_dir = os.getcwd()
    
    analysis_dir = None
    if "base_dir" in analysis_cfg:
        analysis_dir = os.path.join(analysis_cfg["base_dir"], run_id)
    if analysis_dir is None or not os.path.exists(analysis_dir):
        analysis_dir = tempfile.mkdtemp()
        
    dirs = {"work": os.path.normpath(analysis_dir),
            "flowcell": os.path.normpath(archive_dir)}
    assert os.path.exists(dirs["flowcell"]), \
    "The flowcell directory, {}, could not be found".format(dirs["flowcell"])
    assert os.path.exists(dirs["work"]), \
    "The work directory, {}, could not be found".format(dirs["work"])

    if run_info_file is None:
        run_info_file = os.path.join(dirs["flowcell"], "run_info.yaml")
        
        if not os.path.exists(run_info_file):
            # Locate the samplesheet and convert to yaml
            samplesheet = _find_samplesheet(dirs["flowcell"])
            assert samplesheet, \
            "Could not locate samplesheet in {}, aborting..".format(dirs["flowcell"])
            fh, run_info_file = tempfile.mkstemp()
            os.close(fh)
            run_info_file = ssheet.csv2yaml(samplesheet,run_info_file)
            
    assert os.path.exists(run_info_file), \
    "The run info configuration file, {}, could not be found".format(run_info_file)
    
    fc_name, fc_date = fc.get_flowcell_info(dirs["flowcell"])
    # If we have no bc_metrics files in the workdir, we may be looking at a Casava run.
    # In that case, attempt to parse the Demultiplex_Stats.htm file and create bc_metrics files
    metric_files = glob.glob(os.path.join(dirs["work"], "*_barcode", "*bc[_.]metrics")) + glob.glob(os.path.join(dirs["work"], "*bc[_.]metrics"))
    if len(metric_files) == 0:
        casava_report = _find_casava_report(dirs["flowcell"])
        assert len(casava_report) > 0, \
        "Could not locate CASAVA demultiplex report in {}, aborting..".format(dirs["flowcell"])
        metric_files = _casava_report_to_metrics(run_info_file, casava_report, dirs)

    assert len(metric_files) > 0, \
    "Could not locate or create required metric files, aborting.."
    
    print("A report will be created on Google Docs based on the demultiplexed data in {}".format(dirs["work"]))
    print("The configuration file is {0} and the run info file is {1}".format(config_file, run_info_file))
    print("The run was started on {0} and has flowcell id {1}".format(fc_date, fc_name))

    if not dryrun:
        create_report_on_gdocs(fc_date, fc_name, run_info_file, dirs, config)
    else:
        print("DRY-RUN: nothing uploaded")


def _find_casava_report(fc_dir):
    """Locate the CASAVA demultiplex report under the root directory
    of an illumina flowcell output directory
    """
    
    fc_name, _ = fc.get_flowcell_info(fc_dir)
    
    casava_report_glob = os.path.join(fc_dir,"Unaligned*","Basecall_Stats_*{}".format(fc_name[1:]),"Demultiplex_Stats.htm")
    return glob.glob(casava_report_glob)
    
def _find_samplesheet(fc_dir):
    """Locate the samplesheet in the root directory
    of an illumina flowcell output directory
    """
    
    fc_name, _ = fc.get_flowcell_info(fc_dir)
    
    for name in (fc_name, fc_name[1:], "SampleSheet"):
        ssheet = os.path.join(fc_dir, "{}.csv".format(name))
        if os.path.exists(ssheet):
            return ssheet 
    return None
    
def _casava_report_to_metrics(run_info_file, casava_report, dirs):
    """Convert the supplied CASAVA demultiplex report into bcbb-style
    metric files, based on the configuration in the run_info_file.
    Metric files are written to the workdir
    """
    
    metric_files = []
    metrics = defaultdict(dict)
    for report in casava_report:
        for lane, data in dmx._parse_demultiplex_stats_htm(report).items():
            for sequence, metric in data.items():
                # Assert that we are not overwriting a previously parsed metric
                assert not (lane in metrics and sequence in metrics[lane]), \
                "Conflicting demultiplex metrics found for lane {} and index {}. " \
                "This means that there are multiple demultiplex results for the same sample. " \
                "Please review and rectify before proceeding!".format(lane,sequence)
                metrics[lane][sequence] = metric
                
    with open(run_info_file) as fh:
        info = yaml.load(fh)

        fc_name, fc_date = fc.get_flowcell_info(dirs["flowcell"])
        for item in info:
            metrics_file = "{}_{}_{}.bc_metrics".format(item["lane"], fc_date, fc_name)
            multiplex = item.get("multiplex", [])
            for plex in multiplex:
                plex["lane"] = item["lane"]

            dmx._write_demultiplex_metrics(multiplex, metrics, os.path.join(dirs["work"], metrics_file))
            metric_files.append(metrics_file)
    return metric_files
        
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--run-info-file", dest="run_info_file", default=None)
    parser.add_option("-n", "--dry-run", dest="dryrun", action="store_true", default=False)
    options, args = parser.parse_args()

    run_id = None
    if len(args) == 2:
        run_id = args[0]
        config_file = args[1]
    else:
        print(__doc__)
        sys.exit()

    main(run_id, config_file,
         options.run_info_file, options.dryrun)
