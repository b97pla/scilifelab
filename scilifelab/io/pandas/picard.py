"""pm picard lib"""
import os
import re
import pandas as pd
from scilifelab.io import index_containing_substring
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

METRICS_TYPES=['align', 'hs', 'dup', 'insert']

def _raw(x):
    return (x, None)

def _convert_input(x):
    if re.match("^[0-9]+$", x):
        return int(x)
    elif re.match("^[0-9,.]+$", x):
        return float(x.replace(",", "."))
    else:
        return str(x)

def _read_picard_metrics(f):
    if not os.path.exists(f):
        LOG.warn("IO failure: no such file {}".format(f))
        return (None, None)
    with open(f) as fh:
        data = fh.readlines()
        # Find histogram line
        i_hist = index_containing_substring(data, "## HISTOGRAM")
        if i_hist == -1:
            i = len(data)
        else:
            i = i_hist
        tmp = [[_convert_input(y) for y in x.rstrip("\n").split("\t")] for x in data[0:i] if not re.match("^[ #\n]", x)]
        metrics = pd.DataFrame(tmp[1:], columns=tmp[0])
        if i_hist == -1:
            return (metrics, None)
        tmp = [[_convert_input(y) for y in x.rstrip("\n").split("\t")] for x in data[i_hist:len(data)] if not re.match("^[ #\n]", x)]
        hist = pd.DataFrame(tmp[1:], columns=tmp[0])
    return (metrics, hist)

# For now: extension maps to tuple (label, description). Label should
# be reused for analysis definitions
EXTENSIONS={'.align_metrics':('align', 'alignment', _read_picard_metrics),
            '.hs_metrics':('hs', 'hybrid selection', _read_picard_metrics),
            '.dup_metrics':('dup', 'duplication metrics', _read_picard_metrics),
            '.insert_metrics':('insert', 'insert size', _read_picard_metrics),
            '.eval_metrics':('eval', 'snp evaluation', _raw)
            }


def read_metrics(f):
    """Read metrics"""
    (_, metrics_type) = os.path.splitext(f)
    d = EXTENSIONS[metrics_type][2](f)
    return d
