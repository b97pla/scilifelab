"""rst module for generating reports."""
import os
import texttable as tt
from datetime import datetime
from scilifelab.log import minimal_logger
from mako.template import Template
from mako.exceptions import RichTraceback

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
    
def make_rst_sample_table(data):
    """Format sample table"""
    if data is None:
        return ""
    else:
        tab_tt = tt.Texttable()
        tab_tt.add_rows(data)
        return tab_tt.draw()

def make_rest_note(outfile, sample_table=None, outdir="rst", report="sample_report", **kw):
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
               'sample_table':make_rst_sample_table(sample_table),
               'stylefile': os.path.join(TEMPLATEPATH, "rst", "scilife.txt"),
               'sll_logo_small':sll_logo_small,
               })
    if not os.path.exists(sphinx_path):
        os.makedirs(sphinx_path)
    # add makefile if not present
    _install_makefile(sphinx_path, **kw)
    # Write report note
    with open(os.path.join(outdir, os.path.basename(outfile)), "w") as fh:
        fh.write(_render(report_templates[report], **kw))
    
