"""rst module for generating reports."""
import os
import re
import texttable as tt
from datetime import datetime
from scilifelab.log import minimal_logger
from mako.template import Template
from mako.exceptions import RichTraceback
from cStringIO import StringIO

# Set minimal logger
LOG = minimal_logger(__name__)

# Now
now = datetime.now()

# Set sll_logo file
FILEPATH=os.path.dirname(os.path.realpath(__file__))
sll_logo = os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo.gif")
sll_logo_small = os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo_small.gif")

TEMPLATEPATH=os.path.join(FILEPATH, os.pardir, "data", "templates")

report_templates = {'project_report':Template(filename=os.path.join(TEMPLATEPATH, "project_report.mako")),
                    'sample_report':Template(filename=os.path.join(TEMPLATEPATH, "sample_report.mako")),
                    'bp_seqcap':Template(filename=os.path.join(TEMPLATEPATH, "bp_seqcap.mako")),
                    }
rst_templates = {
                 'make':Template(filename=os.path.join(TEMPLATEPATH, "rst", "Makefile.mako"))
                 }


def _render(tpl, **kw):
    """Render a mako template, catching exceptions if present.

    :param tpl: mako template 
    :param kw: keyword arguments for formatting template
    
    """
    res = ""
    try:
        res = tpl.render(**kw)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            LOG.info("File {}, line {}, in {}".format(filename, lineno, function))
            LOG.info(line, "\n")
        LOG.info("{}: {}".format(str(traceback.error.__class__.__name__), traceback.error))
    return res

def _install_makefile(path, **kw):
    """Install rst Makefile to path.

    :param path: path
    :param kw: keyword arguments for formatting
    """
    if not os.path.exists(os.path.join(path, "Makefile")):
        with open(os.path.join(path, "Makefile"), "w") as fh:
            fh.write(_render(rst_templates['make'], **kw))

def indent_texttable_for_rst(ttab, indent=4, add_spacing=True):
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
    
def col_widths(data):
    if not data or len(data) == 0:
        return []
    wd = [0 for i in xrange(len(data[0]))]
    for row in data:
        for i in xrange(len(row)):
            wd[i] = max(len(str(row[i])),wd[i])
    return wd

def make_rst_sample_name_table(data):
    """Format sample name table"""
    if data is None or len(data) == 0:
        return ""
    else:
        tab_tt = tt.Texttable()
        tab_tt.add_rows(data)
        tab_tt.set_cols_width(col_widths(data))
        return indent_texttable_for_rst(tab_tt)

def make_rst_sample_yield_table(data):
    """Format sample name table"""
    if data is None or len(data) == 0:
        return ""
    else:
        tab_tt = tt.Texttable()
        tab_tt.set_precision(2)
        tab_tt.add_rows(data)
        tab_tt.set_cols_width(col_widths(data))
        return indent_texttable_for_rst(tab_tt)

def make_rst_sample_quality_table(data):
    """Format sample name table"""
    if data is None or len(data) == 0:
        return ""
    else:
        tab_tt = tt.Texttable()
        tab_tt.set_precision(2)
        tab_tt.add_rows(data)
        tab_tt.set_cols_width(col_widths(data))
        return indent_texttable_for_rst(tab_tt)

def make_rest_note(outfile, tables={}, outdir="rst", report="sample_report", **kw):
    """Make reSt-formatted note.

    :param outfile: outfile name
    :param sample_table: list of list sample table representation
    :param outdir: output directory
    :param report: default report template to use
    :param kw: keyword arguments
    
    """
    rst_path = os.path.join(os.path.dirname(outfile), outdir)
    kw.update({'date':now.strftime("%Y-%m-%d"),
               'year': now.year,
               'project_lc':kw.get('project_name'),
               'author':'N/A',
               'description':'N/A',
               'sample_name_table':make_rst_sample_name_table(tables.get('name',None)),
               'sample_yield_table':make_rst_sample_yield_table(tables.get('yield',None)),
               'sample_quality_table':make_rst_sample_quality_table(tables.get('quality',None)),
               'stylefile': os.path.join(TEMPLATEPATH, "rst", "scilife.txt"),
               'sll_logo_small':sll_logo_small,
               })
    if not os.path.exists(rst_path):
        os.makedirs(rst_path)
    # add makefile if not present
    _install_makefile(rst_path, **kw)
    # Write report note
    with open(os.path.join(outdir, os.path.basename(outfile)), "w") as fh:
        fh.write(_render(report_templates[report], **kw))
    
def make_sample_rest_notes(concat_outfile, s_param_list, outdir="rst"):
    """Make reST-formatted sample note and concatenated ditto

    :param outdir: output directory
    :param concat_outfile: concatenated outfile
    :param s_param_list: list of samples parameter dictionaries
    """
    concatenated_rst = StringIO()
    for s_param in s_param_list:
        if s_param["outfile"].endswith(".pdf"):
            s_param["outfile"] = s_param["outfile"].replace(".pdf", ".rst")
        rst_path = os.path.join(os.path.dirname(s_param["outfile"]), outdir)
        s_param.update({'date':"{:%B %d, %Y}".format(datetime.now()),
                        'year': now.year,
                        'project_lc':s_param.get('project_name'),
                        'author':'N/A',
                        'description':'N/A',
                        'stylefile': os.path.join(TEMPLATEPATH, "rst", "scilife.txt"),
                        'sll_logo_small':sll_logo_small,
                        })
        if not os.path.exists(rst_path):
            os.makedirs(rst_path)
        # add makefile if not present
        _install_makefile(rst_path, **s_param)
        # Write report note
        rst_out = _render(report_templates["sample_report"], **s_param)
        with open(os.path.join(outdir, os.path.basename(s_param["outfile"])), "w") as fh:
            fh.write(rst_out)
        concatenated_rst.write(rst_out)
        concatenated_rst.write(".. raw:: pdf\n\n   PageBreak\n\n")

    # Write concatenated outfile
    with open(os.path.join(outdir, concat_outfile), "w") as fh:
        fh.write(concatenated_rst.getvalue())
    
