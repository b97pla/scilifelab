"""
rl.py

Format reportlab documents
"""

import sys
import os
from datetime import datetime

from collections import OrderedDict
from mako.template import Template

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Image
from reportlab.rl_config import defaultPageSize

# Set sll_logo file
FILEPATH=os.path.dirname(os.path.realpath(__file__))
sll_logo = os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo.gif")


styles = getSampleStyleSheet()
p = styles['Normal']
h1 = styles['Heading1']
h2 = styles['Heading2']
h3 = styles['Heading3']
h4 = styles['Heading4']

## FIXME: should mako templates go to data/templates?
def sample_note_paragraphs():
    """Get paragraphs for sample notes."""
    paragraphs = OrderedDict()
    paragraphs["Project name"] = dict(style=h3, 
                                      tpl=Template("${project_name} (${customer_reference})"))
    
    paragraphs["UPPNEX project id"] = dict(style=h3, 
                                           tpl=Template("${uppnex_project_id}"))
    
    paragraphs["Flow cell id"] = dict(style=h3, 
                                      tpl=Template("${FC_id}"))
    
    paragraphs["Sequence data directory"] = dict(style=h3, 
                                                 tpl=Template("/proj/${uppnex_project_id}/INBOX/${project_name}/${scilifelab_name}/${start_date}_${FC_id}"))
    
    paragraphs["Sample"] = dict(style=h3,
                                tpl=Template("""${scilifelab_name} / ${customer_name}.
Ordered amount: ${ordered_amount} million paired reads."""))
    
    paragraphs["Method"] = dict(style=h3,
                                tpl = Template("""Clustered on cBot and sequenced on HiSeq 2000
according to manufacturer's instructions. Base
conversion using OLB v1.9, demultiplexed and
converted to fastq using CASAVA v1.8. The quality scale
is Sanger / phred33 / Illumina 1.8+."""))
    
    paragraphs["Results"] = dict(style=h3,
                                 tpl = Template("""${rounded_read_count} million reads in lane with PhiX
error rate ${phix_error_rate}%. Average quality score
${avg_quality_score}."""))
    
    paragraphs["Comments"] = dict(style=h3,
                                  tpl = Template("${success}"))
                                  
    return paragraphs

def formatted_page(canvas, doc):
    """Page format for a document, which adds headers and footers and the like
    which should be on every page of the delivery note.
    """
    canvas.saveState()
    canvas.drawImage(str(sll_logo), 2 * cm, defaultPageSize[1] - 2 * cm, 4 * cm, 1.25 * cm)
    canvas.restoreState()

def make_note(outfile, headers, paragraphs, **kw):
    """Builds a pdf file named outfile based on headers and
    paragraphs, formatted according to parameters in kw.

    :param outfile: outfile name
    :param headers: <OrderedDict> of headers
    :param paragraphs: <OrderedDict> of paragraphs
    :param kw: keyword arguments for formatting
    """
    story = []
    [story.append(Paragraph(x, headers[x])) for x in headers.keys()]

    for headline, paragraph in paragraphs.items():
        story.append(Paragraph(headline, paragraph.get("style", h3)))
        if isinstance(paragraph.get("tpl"), dict):
            for sub_headline, sub_paragraph in paragraph.items():
                story.append(Paragraph(sub_headline, paragraph.get("style", h4)))
                story.append(Paragraph(sub_paragraph.get("tpl").render(**kw),  p))
        else:
            story.append(Paragraph(str(paragraph.get("tpl").render(**kw)), p))

    doc = SimpleDocTemplate(outfile)
    doc.build(story, onFirstPage=formatted_page, onLaterPages=formatted_page)
    return doc

def make_example_note(outfile):
    """Make a note with some simple nonsensical data. Looking at this function
    and running it to make a PDF should give an idea about the structure of the
    script.
    """
    headers = OrderedDict()
    headers["Raw data delivery note"] = h1
    headers["SciLifeLab Stockholm"]   = h2
    headers["{:%B %d, %Y}".format(datetime.now())] = h2
    paragraphs = sample_note_paragraphs()
    kw = {
        "project_name": "A_test",
        "customer_reference": "Some_test",
        "uppnex_project_id": "b2013444",
        "ordered_amount":"23",
        "start_date": "000101",
        "FC_id": "SN001_001_AABCD99XX",
        "scilifelab_name": "Test sample",
        "customer_name": "That sample for a test",
        "rounded_read_count": "1",
        "phix_error_rate": "1",
        "avg_quality_score": "1",
        "success": "How should I know if it was successful or not?",
        }

    make_note(outfile, headers, paragraphs,**kw)

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

