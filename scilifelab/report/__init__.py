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

styles = getSampleStyleSheet()
p = styles['Normal']
h1 = styles['Heading1']
h2 = styles['Heading2']
h3 = styles['Heading3']
h4 = styles['Heading4']

paragraphs = OrderedDict()

paragraphs["abstract"] = Template("""
The ${enrichment_kit} target enrichment kit was used to
prepare and sequence DNA enriched for exon sequence. Sequence data was
mapped to a reference, and variants were called using a publically
available software pipeline.
""")

def halo_report():
    kw = {enrichment_kit= "Agilent SureSelect All Exon 50mb"}
    print kw
    return paragraphs["abstract"].render(**kw)
