"""
ls library
"""

import sys
import csv

def runinfo_to_tab(runinfo_yaml):
    return _yaml_to_tab(runinfo_yaml)

def runinfo_dump(runinfo_tab, fh=sys.stdout):
    w=csv.writer(fh, delimiter="\t")
    w.writerows(runinfo_tab)

def runinfo_projects(runinfo_tab):
    return list(set(_column(runinfo_tab, "sample_prj")))

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

def _column(matrix, label):
    i = matrix[0].index(label)
    return [row[i] for row in matrix[1:]]


