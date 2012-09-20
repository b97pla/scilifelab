"""
runinfo file utilities
"""

import os
import sys
import yaml
import json
import csv
import glob
from cStringIO import StringIO

## FIXME: rewrite this hack
def load_runinfo(config, flowcell):
    if os.path.exists(os.path.join(config.get("archive", "root"), flowcell, "run_info.yaml")):
        runinfo_tab = get_runinfo(os.path.join(config.get("archive", "root"), flowcell, "run_info.yaml"))
    elif os.path.exists(os.path.join(config.get("analysis", "root"), flowcell, "run_info.yaml")):
        runinfo_tab = get_runinfo(os.path.join(config.get("analysis", "root"), flowcell, "run_info.yaml"))
    else:
        self.log.warn("No run information available")
        return None
    return runinfo_tab

def get_runinfo(path, tab=True):
    """Get runinfo contents as tab for a given path."""
    with open(path) as fh:
        runinfo_yaml = yaml.load(fh)
    if tab:
        return _runinfo_to_tab(runinfo_yaml)
    else:
        return runinfo_yaml

def _runinfo_to_tab(runinfo_yaml):
    """Convert yaml to tabular format"""
    return _yaml_to_tab(runinfo_yaml)

def _yaml_to_tab(runinfo_yaml):
    """Convert yaml to tab"""
    out = []
    lanelabels = ["lane", "description", "flowcell_id", "genome_build", "analysis"]
    mplabels = ["barcode_type", "barcode_id", "sample_prj", "name", "sequence"]
    header = lanelabels + mplabels
    out.append(header)
    for info in runinfo_yaml:
        laneinfo = [info.get(x, None) for x in lanelabels]
        for mp in info.get("multiplex", None):
            mpinfo = [mp.get(x, None) for x in mplabels]
            line = laneinfo + mpinfo
            out.append(line)
    return out

def _get_rows(matrix, i):
    return [row for row in matrix[i:]]

def _column(matrix, label):
    i = matrix[0].index(label)
    return [row[i] for row in matrix[1:]]

def _convert_barcode_id_to_name(multiplex, fc_name, fq):
    bcid2name = dict([(str(mp.get_barcode_id()), mp.get_full_name()) for mp in multiplex])
    if options.no_full_names:
        bcid2name = dict([(str(mp.get_barcode_id()), get_sample_name(mp.get_barcode_name())) for mp in multiplex])
    bcid = re.search("_(\d+)_(\d+)_fastq.txt", fq)
    from_str = "%s_%s_fastq.txt" % (bcid.group(1), bcid.group(2))
    to_str   = "%s_%s.fastq" % (bcid2name[bcid.group(1)], bcid.group(2))
    return fq.replace(from_str, to_str)


def dump_runinfo(runinfo, tab=True):
    """Dump runinfo tabular information to output """
    fh = StringIO()
    if tab:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(runinfo)
    else:
        fh.write(json.dumps(runinfo))
    return fh.getvalue()

def runinfo_projects(runinfo_tab):
    """List runinfo projects"""
    return list(set(_column(runinfo_tab, "sample_prj")))

def runinfo_lanes(runinfo_tab):
    """List runinfo lanes"""
    return list(set(_column(runinfo_tab, "lane")))

def runinfo_barcodes(runinfo_tab, lane):
    """List runinfo barcodes for a lane"""
    info = subset_runinfo(runinfo_tab, "lane", lane)
    return list(set(_column(info, "barcode_id")))

def subset_runinfo(runinfo_tab, column, query):
    """Get subset of runinfo by subsetting a header"""
    vals = list(_column(runinfo_tab, column))
    i = [0] + [j + 1 for j in range(0, len(vals)) if vals[j]==query]
    return [runinfo_tab[j] for j in i ]

def find_files(runinfo_tab, path):
    glob_str = [os.path.join(path, "{}_*_?{}_barcode/{}_[0-9][0-9][0-9][0-9][0-9][0-9]_?{}_nophix_{}*".format(x[0], x[2], x[0], x[2], x[6])) for x in runinfo_tab[1:]]
    print glob_str
    res = [glob.glob(x) for x in glob_str]
    return res


# def list_runinfo(self):
#     """List runinfo for a given flowcell"""
#     if self.pargs.flowcell is None or self.pargs.flowcell == "default":
#         return "Please provide flowcell id"
#     assert self.config.get("archive", "root"), "archive directory not defined"
#     f = os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml")
#     self.log.info("Opening file %s" %f)
#     with open(f) as fh:
#         runinfo_yaml = yaml.load(fh)
#     runinfo_tab = runinfo_to_tab(runinfo_yaml)
#     if self.pargs.tab:
#         return runinfo_dump(runinfo_tab)
#     elif self.pargs.list_projects:
#         return "available projects for flowcell %s:\n\t" %self.pargs.flowcell + "\n\t".join(runinfo_projects(runinfo_tab))
#     else:
#         return runinfo_yaml

