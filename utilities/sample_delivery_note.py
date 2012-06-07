from datetime import datetime
from collections import OrderedDict

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate
# from reportlab.platypus import Image
# from reportlab.platypus import Spacer

from reportlab.rl_config import defaultPageSize

styles = getSampleStyleSheet()
p = styles['Normal']
h1 = styles['Heading1']
h2 = styles['Heading2']
h3 = styles['Heading3']
h4 = styles['Heading4']

paragraphs = OrderedDict()
paragraphs["Project name"] = "{project_name} ({customer_reference})"
paragraphs["UPPNEX project id"] = "{uppnex_project_id}"
paragraphs["Sequence data directory"] = \
"/proj/{uppnex_project_id}/INBOX/{project_name}/{start_date}_{FC_id}"
paragraphs["Sample"] = "{scilifelab_name} / {customer_name}\n\n" \
                       "Ordered amount in millions of read (pairs)."
paragraphs["Method"] = "Clustered on cBot and sequenced on HiSeq 2000 " \
                       "according to manufacturer's instructions. Base " \
                       "conversion using OLB v1.9, demultiplexed and " \
                       "converted to fastq using CASAVA v1.8."
paragraphs["Results"] = "{rounded_read_count} million reads in lane with PhiX " \
                        "error rate {phix_error_rate}%. Average quality score " \
                        "{avg_quality_score}."
paragraphs["Comments"] = "Successful run/unsuccessful run, enough data/not " \
                         "enough data"
paragraphs["Information"] = OrderedDict()
paragraphs["Information"]["Acknowledgement"] = \
"Please notify us when you publish using data produced at Science For Life " \
"Laboratory (SciLifeLab) Stockholm. To acknowledge SciLifeLab Stockholm in your " \
'article, you can use a sentence like "The authors would like to acknowledge ' \
"support from Science for Life Laboratory, the national infrastructure SNISS, " \
"and Uppmax for providing assistance in massively parallel sequencing and " \
'computational infrastructure."'
paragraphs["Information"]["Naming conventions"] = \
"The data is delivered in fastq format using Illumina 1.8 quality scores. " \
"There will be one file for the forward reads and one file for the reverse " \
"reads. More information on our naming conventions can be found here."
paragraphs["Information"]["Data access at UPPMAX"] = \
"Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next " \
"Generation sequence Cluster & Storage, www.uppmax.uu.se), from which the " \
"user can access it. If you have problems to access your data, please contact " \
"SciLifeLab genomics_support@scilifelab.se. If you have questions regarding " \
"UPPNEX, please contact support@uppmax.uu.se. Information on how to access your " \
"data can be found here."


def formatted_page(canvas, doc):
    canvas.saveState()
    canvas.drawImage("sll_logo.gif", 2 * cm, defaultPageSize[1] - 2 * cm, 4 * cm, 1.25 * cm)
    canvas.restoreState()


def make_note(parameters):
    story = []
    story.append(Paragraph("Raw data delivery note", h1))
    story.append(Paragraph("SciLifeLab Stockholm", h2))
    story.append(Paragraph("{:%B %d, %Y}".format(datetime.now()), h2))

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


"""
story.append(Paragraph("Project name", h3))
story.append(Paragraph("A_test (Some_test)", p))
story.append(Paragraph("UPPNEX project id", h3))
story.append(Paragraph("b2013444", p))
story.append(Paragraph("Sequence data directory", h3))
story.append(Paragraph(directory, p))
story.append(Paragraph("Sample", h3))
story.append(Paragraph(sample, p))
story.append(Paragraph("Method", h3))
story.append(Paragraph(method, p))
story.append(Paragraph("Results", h3))
story.append(Paragraph(results, p))
story.append(Spacer(0, 4 * cm * 2 + 1 * cm))  # Graphs should go here
story.append(Paragraph("Comments", h3))
story.append(Paragraph(comments, p))
story.append(Paragraph("Information", h3))
story.append(Paragraph("Acknowledgement", h4))
story.append(Paragraph(acknowledgement, p))
story.append(Paragraph("Naming conventions", h4))
story.append(Paragraph(naming, p))
story.append(Paragraph("Data access at UPPMAX", h4))
story.append(Paragraph(data_access, p))
doc = SimpleDocTemplate("mydoc.pdf")
doc.build(story, onFirstPage=formatted_page, onLaterPages=formatted_page)
"""


def make_example_note():
    """Make a note with some simple nonsensical data.
    """
    parameters = {
    "project_name": "A_test",
    "customer_reference": "Some_test",
    "uppnex_project_id": "b2013444",
    "start_date": "000101",
    "FC_id": "SN001_001_AABCD99XX",
    "scilifelab_name": "Test sample",
    "customer_name": "That sample for a test",
    "rounded_read_count": "1",
    "phix_error_rate": "1",
    "avg_quality_score": "1",
    }

    make_note(parameters)
