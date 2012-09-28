"""
Reporting utilities
"""

import sys
from mako.template import Template
from collections import OrderedDict

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
from reportlab.rl_config import defaultPageSize

FILEPATH=os.path.dirname(os.path.realpath(__file__))

styles = getSampleStyleSheet()
p = styles['Normal']
h1 = styles['Heading1']
h2 = styles['Heading2']
h3 = styles['Heading3']
h4 = styles['Heading4']

def formatted_page(canvas, doc):
    """Page format for a document, which adds headers and footers and the like
    which should be on every page of the delivery note.
    """
    canvas.saveState()
    canvas.drawImage(os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo.gif"), 2 * cm, defaultPageSize[1] - 2 * cm, 4 * cm, 1.25 * cm)
    canvas.restoreState()

def make_note(**header, **paragraphs, **kw):
    """Builds a pdf note based on the passed dictionary of parameters.
    For the structure of the parameters, see the code of the function
    make_example_note.
    """
    story = []
    [story.append(Paragraph(x, header[x]))]

    for headline, paragraph in paragraphs.items():
        story.append(Paragraph(headline, h3))
        if isinstance(paragraph, dict):
            for sub_headline, sub_paragraph in paragraph.items():
                story.append(Paragraph(sub_headline, h4))
                story.append(Paragraph(sub_paragraph, p))
        else:
            story.append(Paragraph(paragraph.format(**parameters), p))

    doc = SimpleDocTemplate("{scilifelab_name}_note.pdf".format(**parameters))
    doc.build(story, onFirstPage=formatted_page, onLaterPages=formatted_page)

def calc_avg_qv(counts):
    qv = 2
    qsum = 0
    seqsum = 0
    for v in counts:
        qsum += qv * float(v)
        seqsum += float(v)
        qv += 1
    avg_q = qsum / seqsum
    return round(avg_q, 1)

