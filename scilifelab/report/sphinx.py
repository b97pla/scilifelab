"""Sphinx module for generating reports."""
import os
from datetime import datetime
from scilifelab.log import minimal_logger
from mako.template import Template

# Set minimal logger
LOG = minimal_logger(__name__)

# Now
now = datetime.now()

# Set sll_logo file
FILEPATH=os.path.dirname(os.path.realpath(__file__))
sll_logo = os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo.gif")

TEMPLATEPATH=os.path.join(FILEPATH, os.pardir, "data", "templates")

report_templates = {'project_report':Template(filename=os.path.join(TEMPLATEPATH, "project_report.mako")),
                    'sample_report':Template(filename=os.path.join(TEMPLATEPATH, "sample_report.mako")),
                    'bp_seqcap':Template(filename=os.path.join(TEMPLATEPATH, "bp_seqcap.mako")),
                    }
sphinx_templates = {'conf':Template(filename=os.path.join(TEMPLATEPATH, "sphinx", "conf.mako")),
                    'make':Template(filename=os.path.join(TEMPLATEPATH, "sphinx", "Makefile.mako"))
                    }


def _install_sphinx_files(path, **kw):
    """Install sphinx conf.py and Makefile to path.

    :param path: path
    :param kw: keyword arguments for formatting
    """
    if not os.path.exists(os.path.join(path, "Makefile")):
        with open(os.path.join(path, "Makefile"), "w") as fh:
            fh.write(sphinx_templates['make'].render(**kw))
    if not os.path.exists(os.path.join(path, "conf.py")):
        with open(os.path.join(path, "conf.py"), "w") as fh:
            fh.write(sphinx_templates['conf'].render(**kw))

def make_sphinx_sample_table(data):
    """Format sample table"""
    pass

def make_rest_note(outfile, outdir="sphinx", report="sample_report", **kw):
    """Make reSt-formatted note."""
    sphinx_path = os.path.join(os.path.dirname(outfile), outdir)
    kw.update({'date':now.strftime("%Y-%m-%d"),
               'year': now.year,
               'project_lc':kw.get('project_name'),
               'author':'N/A',
               'description':'N/A',
               })
    if not os.path.exists(sphinx_path):
        os.makedirs(sphinx_path)
        _install_sphinx_files(sphinx_path, **kw)
    with open(os.path.join(outdir, os.path.basename(outfile)), "w") as fh:
        fh.write(report_templates[report].render(**kw))
    
