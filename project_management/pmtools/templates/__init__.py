"""
Templates configuration
"""

import os
import sys
from mako.template import Template
from mako.lookup import TemplateLookup

TEMPLATE_ROOT = os.path.join(os.path.abspath(os.path.dirname(__file__)), "tpl")
MAKE_TEMPLATE_DIR = os.path.join(TEMPLATE_ROOT, "make")

def get_make_templates():
    tmpl = {'Makefile' : Template(filename=os.path.join(MAKE_TEMPLATE_DIR, 'Makefile'))}
    return tmpl
