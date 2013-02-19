"""Module best_practice - code for generating best practice reports and notes"""
import os
import re
import yaml
import pandas as pd
from cStringIO import StringIO
from scilifelab.report.rst import make_rest_note
from texttable import Texttable
from itertools import izip
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

BEST_PRACTICE_NOTES=["seqcap"]

SEQCAP_TABLE_COLUMNS = ["Sample", "Total", "Aligned", "Pair duplicates", "Insert size", "On target", "Mean coverage", "10X coverage", "0X coverage", "Variations", "In dbSNP", "Ts/Tv (all)", "Ts/Tv (dbSNP)", "Ts/Tv (novel)"]

SEQCAP_KITS={
    'agilent_v4':'Agilent SureSelect XT All Exon V4',
    'agilent_v5':'Agilent SureSelect Human All Exon V5',
    'agilent_v5_utr':'Agilent SureSelect Human All Exon V5 UTRs',
    'custom':'Custom',
    }

parameters = {
    'projectsummarytable' : None,
    'projecttableref' : None,
    }

def _dataframe_to_texttable(df, align=None):
    """Convert data frame to texttable. Sets column widths to the
    widest entry in each column."""
    ttab = Texttable()
    ttab.set_precision(1)
    h = [[x for x in df]]
    h.extend([x for x in df.to_records(index=False)])
    if align:
        colWidths = [max(len(x), len(".. class:: {}".format(y))) for x,y in izip(df.columns, align)]
    else:
        colWidths = [len(x) for x in df.columns]
    for row in h:
        for i in range(0, len(row)):
            if type(row[i]) == str:
                colWidths[i] = max([len(str(x)) for x in row[i].split("\n")] + [colWidths[i]])
            colWidths[i] = max(len(str(row[i])), colWidths[i])
    table_data = []
    if align:
        for row in h:
            table_row = []
            i = 0
            for col, aln in izip(row, align):
                table_row.append(".. class:: {}".format(aln) + " " * colWidths[i] + "{}".format(col))
                i = i + 1
            table_data.append(table_row)
    else:
        table_data = h
    ttab.add_rows(table_data)
    ttab.set_cols_width(colWidths)
    # Note: this does not affect the final pdf output
    ttab.set_cols_align(["r"] * len(colWidths))
    return ttab    

def _indent_texttable_for_rst(ttab, indent=4, add_spacing=True):
    """Texttable needs to be indented for rst.

    :param ttab: texttable object
    :param indent: indentation (should be 4 *spaces* for rst documents)
    :param add_spacing_row: add additional empty row below class directives

    :returns: reformatted texttable object as string
    """
    output = ttab.draw()
    new_output = []
    for row in output.split("\n"):
        new_output.append(" " * indent + row)
        if re.search('.. class::', row):
            new_row = [" " if x != "|" else x for x in row]
            new_output.append(" " * indent + "".join(new_row))
    return "\n".join(new_output)


def _split_project_summary_sample_name(samplename):
    """Project summary name consists of description;lane;sequence.
    Description in turn is made up of {sample_prj}_{name}."""
    info = {'description':None, 'Lane':None, 'Sequence':None, 'sample_prj':None, 'Sample':samplename, 'ScilifeName':samplename, 'CustomerName':None}
    if samplename.count(";") == 2:
        info['description'] = samplename.split(";")[0]
        info['Lane'] = samplename.split(";")[1]
        info['Sequence'] = samplename.split(";")[2]
    if info['description']:
        m = re.match("([0-9A-Za-z\.\_]+)_(P[0-9][0-9][0-9]_[0-9A-Za-z\_]+)", info['description'])
        if m.groups():
            info['sample_prj'] = m.groups()[0]
            info['Sample'] = m.groups()[1]
            info['ScilifeName'] = m.groups()[1]
    return info

def _get_seqcap_summary(flist):
    """Gather relevant information for sequence capture."""
    df_list = []
    for run_info in flist:
        prj_summary = os.path.join(os.path.dirname(run_info), "project-summary.csv")
        if not os.path.exists(prj_summary):
            LOG.warn("No project summary file for {}: skipping".format(os.path.basename(run_info)))
            continue
        with open(prj_summary) as fh:
            LOG.debug("Reading file {}".format(prj_summary))
            df_list.append(pd.io.parsers.read_csv(fh, sep=","))
    df = pd.concat(df_list)
    samples_list = [_split_project_summary_sample_name(x) for x in df["Sample"]]
    samples_df = pd.DataFrame([_split_project_summary_sample_name(x) for x in df["Sample"]])
    df["Sample"] = [_split_project_summary_sample_name(x)['Sample'] for x in df["Sample"]]
    df.columns = SEQCAP_TABLE_COLUMNS
    return df.sort(["Sample"]), samples_df.sort(["Sample"])

def _get_software_table(flist):
    df_list = []
    for run_info in flist:
        software_versions = os.path.join(os.path.dirname(run_info), "bcbb_software_versions.txt")
        if not os.path.exists(software_versions):
            LOG.warn("No software versions file for {}: skipping".format(os.path.basename(run_info)))
            continue
        with open(software_versions) as fh:
            LOG.debug("Reading file {}".format(software_versions))
            df_list.append(pd.io.parsers.read_table(fh, sep="\t", skiprows=1, header=None))
    ret_df = df_list[0]
    ret_df.columns = ["Software", "Version"]
    return ret_df.sort(["Software"])

def _get_database_table(flist, post_process=None):
    """Get database versions"""
    if post_process:
        with open(post_process) as fh:
            conf = yaml.load(fh)
    df_list = []
    for run_info in flist:
        with open(run_info) as fh:
            bcbb_conf = yaml.load(fh)
        analyses = list(set([x.get("analysis", None) for x in bcbb_conf.get("details", {})]))
        if len(analyses) > 1:
            LOG.warn("More than one defined analysis: {}\nArbitrarily using the first.".format(analyses)
)
        if not post_process:
            with open(run_info.replace("bcbb-config.yaml", "post_process.yaml")) as fh:
                conf = yaml.load(fh)
        tab = conf.get("custom_algorithms", {}).get(analyses.pop(), {})
        df_list.append(pd.DataFrame({'Database':tab.keys(), 'Version':tab.values()}))

    ret_df = df_list[0]
    ret_df.index = ret_df["Database"]
    i = [x in ["dbsnp", "train_hapmap", "train_indels", "train_1000g_omni"] for x in ret_df.index]
    ret_df = ret_df[i]
    ret_df["Version"] = [os.path.basename(x) for x in ret_df["Version"]]
    return ret_df.sort(["Database"])
    
def _format_num_reads(reads):
    """Format number of reads as k, M, or G"""
    fmt_reads = []
    for x in reads:
        if int(x) < 1e6:
            fmt_reads.append("{:.1f}k".format(int(x)/1e3))
        elif int(x) > 1e9:
            fmt_reads.append("{:.1f}G".format(int(x)/1e9))
        else:
            fmt_reads.append("{:.1f}M".format(int(x)/1e6))
    return fmt_reads
    

def best_practice_note(project_name=None, samples=None, capture_kit="agilent_v4", application="seqcap", flist=[], sample_name_map=None, **kw):
    """Make a best practice application note.

    NB: currently only works for seqcap application.

    :param project_name: project name
    :param samples: samples to work on. Defaults to all samples.
    :param application: chosen application
    """
    param = parameters
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    if application not in BEST_PRACTICE_NOTES:
        LOG.warn("No such application '{}'. Valid choices are: \n\t{}".format(application, "\n\t".join(BEST_PRACTICE_NOTES)))
    if application == "seqcap":
        df, samples_df = _get_seqcap_summary(flist)
        software_df = _get_software_table(flist)
        database_df = _get_database_table(flist, post_process=kw.get("post_process", None))
        if sample_name_map:
            samples_df.CustomerName = [sample_name_map[s]['customer_name'] for s in samples_df.Sample]
        df.Total = _format_num_reads(df.Total)
        ttab = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[1:5]], align=["left", "right", "right", "right", "right"]))
        ttab_target = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[5:9]], align=["left", "right", "right", "right", "right"]))
        ttab_dbsnp = _indent_texttable_for_rst(_dataframe_to_texttable(df[["Sample"] + SEQCAP_TABLE_COLUMNS[9:14]], align=["left", "right", "right", "right", "right", "right"]))
        ttab_samples = _indent_texttable_for_rst(_dataframe_to_texttable(samples_df[["Sample", "CustomerName", "Sequence"]], align=["left", "right", "right"]))
        ttab_software = _indent_texttable_for_rst(_dataframe_to_texttable(software_df, align=["left", "right"]))
        ttab_database = _indent_texttable_for_rst(_dataframe_to_texttable(database_df, align=["left", "right"]))
        param.update({'project_summary':ttab, 'project_target_summary':ttab_target, 'project_dbsnp_summary':ttab_dbsnp, 'table_sample_summary':ttab_samples, 'capturekit':SEQCAP_KITS[capture_kit],
                      'software_versions_table':ttab_software, 'database_versions_table': ttab_database})
        param['project_name'] = project_name if project_name else kw.get("statusdb_project_name", None)
    # Add applications here
    else:
        pass
    # Generic rest call for all templates
    make_rest_note("{}_best_practice.rst".format(kw.get("project", None)), report="bp_seqcap", outdir=kw.get("basedir", os.curdir), **param)
    return output_data
