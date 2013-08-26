"""Reportlab module for generating pdf documents"""
import os
import unicodedata
from datetime import datetime
from pyPdf import PdfFileWriter, PdfFileReader
from collections import OrderedDict
from mako.template import Template

from scilifelab.log import minimal_logger

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image
from reportlab.lib import colors
from reportlab.lib.pagesizes import inch
from reportlab.rl_config import defaultPageSize
from reportlab.pdfbase.pdfmetrics import stringWidth

## Set minimal logger
LOG = minimal_logger(__name__)

# Set sll_logo file
FILEPATH=os.path.dirname(os.path.realpath(__file__))
sll_logo = os.path.join(FILEPATH, os.pardir, "data", "grf", "sll_logo.gif")

styles = getSampleStyleSheet()
p = styles['Normal']
h1 = styles['Heading1']
h2 = styles['Heading2']
h3 = styles['Heading3']
h4 = styles['Heading4']

def sample_note_paragraphs():
    """Get paragraphs for sample notes."""
    paragraphs = OrderedDict()
    paragraphs["Project name"] = dict(style=h3,
                                      tpl=Template("${project_name} ${'({})'.format(customer_reference) if customer_reference not in ['', 'N/A'] else ''}"))

    paragraphs["UPPNEX project id"] = dict(style=h3,
                                           tpl=Template("${uppnex_project_id}"))

    paragraphs["Flow cell id"] = dict(style=h3,
                                      tpl=Template("${FC_id}"))

    paragraphs["Sequence data directory"] = dict(style=h3,
                                                 tpl=Template("/proj/${uppnex_project_id}/INBOX/${project_name}/${scilifelab_name}/${start_date}_${FC_id}"))

    paragraphs["Sample"] = dict(style=h3,
                                tpl=Template("""${scilifelab_name} / ${customer_name}.
Ordered amount: ${ordered_amount} million read${'{}'.format(' pair') if is_paired else ''}s."""))

    paragraphs["Method"] = dict(style=h3,
                                tpl = Template("""Clustered on cBot
and sequenced on ${instrument_version} according to manufacturer's
instructions. Demultiplexing and conversion using ${casava_version}.
The quality scale is Sanger / phred33 / Illumina 1.8+."""))
    paragraphs["Results"] = dict(style=h3,
                                 tpl = Template("""${rounded_read_count} million read${'{}'.format(' pair') if is_paired else ''}s${' in lane with PhiX error rate {}%'.format(phix_error_rate) if phix_error_rate != 'N/A' else ''}.
                                 Average quality score ${avg_quality_score} (${pct_q30_bases}% bases >= Q30)."""))

    paragraphs["Comments"] = dict(style=h3,
                                  tpl = Template("${success}"))
    return paragraphs

def sample_note_headers():
    """Get headers for sample notes."""
    headers = OrderedDict()
    headers["Raw data delivery note"] = h1
    headers["SciLifeLab Stockholm"]   = h2
    headers["{:%B %d, %Y}".format(datetime.now())] = h2
    return headers

def make_sample_table(data,  header_size=11, row_size=9, **kw):
    """Format sample table"""
    default_colWidths=[1.5*inch, 2*inch, 1*inch, 1*inch, 1*inch, .8*inch]

    if not kw.get("colWidths", None):
        colWidths = len(data[0]) * [0*inch]
        for row in data:
            for i in range(0, len(row)):
                width = stringWidth(str(row[i]), "Helvetica", header_size)
                colWidths[i] = min(default_colWidths, max(colWidths[i], 1.2 * width))

    if not kw.get("rowHeights", None):
        rowHeights = len(data)*[0.25*inch]
    t=Table(data, colWidths, rowHeights)

    t.setStyle(TableStyle([
                           ('VALIGN',(0,0),(-1,-1),'BOTTOM'),
                           ('ALIGN',(1,0),(-1,-1),'RIGHT'),
                           ('ALIGN',(0,0),(0,-1),'LEFT'),
                           ('FONTSIZE',(0,0),(-1,0),header_size),
                           ('FONTSIZE',(0,1),(-1,-1), row_size),
                           ('LINEABOVE',(0,1),(-1,1),1,colors.black),
                           ]))
    return t



def project_note_paragraphs():
    """Get paragraphs for project notes."""
    paragraphs = OrderedDict()
    paragraphs["Project name"] = dict(style=h3, tpl=Template("${project_name} ${'({})'.format(customer_reference) if customer_reference not in ['','N/A'] else ''}"))

    paragraphs["UPPNEX project id"] = dict(style=h3, tpl=Template("${uppnex_project_id}"))

    paragraphs["Sequence data directories"] = dict(style=h3, tpl=Template("/proj/${uppnex_project_id}/INBOX/${project_name}/"))

    paragraphs["Samples"] = dict(style=h3, tpl=Template(""))

    paragraphs["Comments"] = dict(style=h3, tpl=Template("${finished}"))

    paragraphs["Information"] = OrderedDict()

    paragraphs["Information"]["Naming conventions"] = dict(
        style=h4,
        tpl=Template("""The data is delivered in fastq format using Illumina 1.8
quality scores. There will be one file for the forward reads and
one file for the reverse reads. More information on our naming
conventions can be found at
http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf."""))

    paragraphs["Information"]["Data access at UPPMAX"] = dict(
        style=h4,
        tpl=Template("""Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next
Generation sequence Cluster & Storage, www.uppmax.uu.se), from which
the user can access it. If you have problems to access your data,
please contact SciLifeLab genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.
Information on how to access your data can be found at
http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf."""))

    paragraphs["Information"]["Acknowledgement"] = dict(
        style=h4,
        tpl=Template("""In publications based on data from the work covered by
this contract, the authors must acknowledge SciLifeLab, NGI and Uppmax: \"The
authors would like to acknowledge support from Science for Life Laboratory,
the National Genomics Infrastructure, NGI, and Uppmax for providing assistance
in massive parallel sequencing and computational infrastructure.\"""" ))
    return paragraphs

def project_note_headers():
    """Get headers for sample notes."""
    headers = OrderedDict()
    headers["Project status note"] = h1
    headers["SciLifeLab Stockholm"]   = h2
    headers["{:%B %d, %Y}".format(datetime.now())] = h2
    return headers

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
    story = [Paragraph(x, headers[x]) for x in headers.keys()]

    #Reportlab does not allow non-ascii characters. This normalizes the string.
    for k, v in kw.iteritems():
        if type(v) == unicode:
            LOG.debug('Param %s' % kw[k])
            kw[k] = unicodedata.normalize('NFKD', v).encode('ascii', 'ignore')

    for headline, paragraph in paragraphs.items():
        story.append(Paragraph(headline, paragraph.get("style", h3)))
        if not paragraph.has_key("tpl"):
            for sub_headline, sub_paragraph in paragraph.items():
                story.append(Paragraph(sub_headline, paragraph.get("style", h4)))
                story.append(Paragraph(sub_paragraph.get("tpl").render(**kw),  p))
        else:
            if isinstance(paragraph.get("tpl"), Template):
                story.append(Paragraph(paragraph.get("tpl").render(**kw), p))
            elif isinstance(paragraph.get("tpl"), Table):
                story.append(Spacer(1, 0.2 * inch))
                story.append(paragraph.get("tpl"))
                story.append(Spacer(1, 0.2 * inch))
            else:
                pass

    doc = SimpleDocTemplate(outfile)
    doc.build(story, onFirstPage=formatted_page, onLaterPages=formatted_page)
    return doc

def concatenate_notes(notes, outfile, numpages=1):
    """Concatenate documents. Warn if numpages in document > numpages.

    :param notes: list of notes
    :param outfile: outfile name
    :param numpages: number of pages to include
    """
    output = PdfFileWriter()
    for n in notes:
        if isinstance(n, SimpleDocTemplate):
            fn = n.filename
        else:
            fn = n
        if not os.path.exists(fn):
            LOG.warn("No such file {}".format(fn))
            continue
        LOG.debug("Adding file {}".format(fn))
        inpdf = PdfFileReader(file(fn, "rb"))
        if inpdf.getNumPages() > numpages:
            LOG.warn("Number of pages in document ({}) larger than number of added pages ({}): not adding complete document".format(inpdf.getNumPages, numpages))
        [output.addPage(inpdf.getPage(i)) for i in range(0, numpages)]
    outputStream = file(outfile, "wb")
    output.write(outputStream)
    outputStream.close()

def make_example_project_note(outfile):
    """Make a note with some simple nonsensical data. Looking at this function
    and running it to make a PDF should give an idea about the structure of the
    script.
    """
    headers = project_note_headers()
    paragraphs = project_note_paragraphs()
    kw = {
    "project_name": "A_test",
    "customer_reference": "Some_test",
    "uppnex_project_id": "b2013444",
    "finished":None,
    }

    LOG.debug("Making example project note with parameters {}".format(kw))
    make_note(outfile, headers, paragraphs, **kw)


def make_example_sample_note(outfile):
    """Make a note with some simple nonsensical data. Looking at this function
    and running it to make a PDF should give an idea about the structure of the
    script.
    """
    headers = sample_note_headers()
    paragraphs = sample_note_paragraphs()
    kw = {
        "project_name": "A_test",
        "customer_reference": "Some_test",
        "is_paired": True,
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
        "instrument":"HiSeq 2000",
        "baseconversion_version":"OLB v1.9",
        "casava_version" : "CASAVA v1.8",
        }

    LOG.debug("Making example sample note with parameters {}".format(kw))
    make_note(outfile, headers, paragraphs, **kw)
