"""Module best_practice - code for generating best practice reports and notes"""
import os
import re
import pandas as pd
from cStringIO import StringIO
from scilifelab.report.rst import make_rest_note
from scilifelab.bcbio.run import find_samples
from texttable import Texttable
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

BEST_PRACTICE_NOTES=["seqcap"]

SEQCAP_TABLE_COLUMNS = ["Sample", "Total", "Aligned", "Pair duplicates", "Insert size", "On target", "Mean coverage", "10X coverage", "0X coverage", "Variations", "In dbSNP", "Ts/Tv (all)", "Ts/Tv (dbSNP)", "Ts/Tv (novel)"]

SEQCAP_KITS={
    'agilent_v4':'Agilent SureSelect XT All Exon V4',
    'agilent_v5':'Agilent SureSelect Human All Exon V5',
    'agilent_v5_utr':'Agilent SureSelect Human All Exon V5 UTRs',
    }

parameters = {
    'projectsummarytable' : None,
    'projecttableref' : None,
    }

def _dataframe_to_texttable(df):
    """Convert data frame to texttable. Sets column widths to the
    widest entry in each column."""
    ttab = Texttable()
    h = [[x for x in df]]
    h.extend([list(x) for x in df.to_records(index=False)])
    ttab.add_rows(h)
    colWidths = [len(x) for x in df.columns]
    for row in df.itertuples(index=False):
        for i in range(0, len(row)):
            colWidths[i] = max(len(str(row[i])), colWidths[i])
    ttab.set_cols_width(colWidths)
    ttab.set_cols_align(["r"] * len(colWidths))
    return ttab    

def _indent_texttable_for_rst(ttab, indent=4):
    """Texttable needs to be indented for rst"""
    return " " * indent + ttab.draw().replace("\n", "\n" + " " * indent)

def _get_seqcap_summary(flist):
    """Gather relevant information for sequence capture."""
    dflist = []
    for run_info in flist:
        prj_summary = os.path.join(os.path.dirname(run_info), "project-summary.csv")
        if not os.path.exists(prj_summary):
            LOG.warn("No project summary file for {}: skipping".format(os.path.basename(run_info)))
            continue
        with open(prj_summary) as fh:
            LOG.debug("Reading file {}".format(prj_summary))
            dflist.append(pd.io.parsers.read_csv(fh, sep=","))
    df = pd.concat(dflist)
    df["Sample"] = [re.sub('.*(P\d+)', '\\1', re.sub(';.*', '', x)) for x in df["Sample"]]
    df.columns = SEQCAP_TABLE_COLUMNS
    return df

def best_practice_note(project_id=None, samples=None, capture_kit="agilent_v4", application="seqcap", flist=[], **kw):
    """Make a best practice application note.

    NB: currently only works for seqcap application.

    :param project_id: project id
    :param samples: samples to work on. Defaults to all samples.
    :param application: chosen application
    """
    param = parameters
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    if application not in BEST_PRACTICE_NOTES:
        LOG.warn("No such application '{}'. Valid choices are: \n\t{}".format(application, "\n\t".join(BEST_PRACTICE_NOTES)))
    if application == "seqcap":
        df = _get_seqcap_summary(flist)
        ttab = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[1:5]]))
        ttab_target = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[6:9]]))
        ttab_dbsnp = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[10:14]]))

    param.update({'project_summary':ttab, 'project_target_summary':ttab_target, 'project_dbsnp_summary':ttab_dbsnp, 'capturekit':SEQCAP_KITS[capture_kit]})
    make_rest_note("{}_best_practice.rst".format(kw.get("project", None)), report="bp_seqcap", outdir=kw.get("basedir", os.curdir), **param)
    return output_data
