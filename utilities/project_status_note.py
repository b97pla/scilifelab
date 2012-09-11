"""Script for generating a PDF project status note for a sample.
The script uses reportlab for PDF, for which a user guide can be found here:
http://www.reportlab.com/software/opensource/rl-toolkit/guide/
"""

import sys 
import couchdb
import optparse 

from datetime import datetime
from collections import OrderedDict
from types import *

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate
# from reportlab.platypus import Image
# from reportlab.platypus import Spacer
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

# This is basically a template of ('headline', 'content') pairs.
# Use dictionary based python string formatting to fill out variables.
paragraphs = OrderedDict()
paragraphs["Project name"] = "{project_name} ({customer_reference})"

paragraphs["UPPNEX project id"] = "{uppnex_project_id}"

paragraphs["Sequence data directories"] = \
"/proj/{uppnex_project_id}/INBOX/{project_name}/"

paragraphs["Samples"] = ""

paragraphs["Comments"] = "{finished}"

paragraphs["Information"] = OrderedDict()

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

paragraphs["Information"]["Acknowledgement"] = \
"Please notify us when you publish using data produced at Science For Life " \
"Laboratory (SciLifeLab) Stockholm. To acknowledge SciLifeLab Stockholm in your " \
'article, you can use a sentence like "The authors would like to acknowledge ' \
"support from Science for Life Laboratory, the national infrastructure SNISS, " \
"and Uppmax for providing assistance in massively parallel sequencing and " \
'computational infrastructure."'


def formatted_page(canvas, doc):
    """Page format for a document, which adds headers and footers and the like
    which should be on every page of the delivery note.
    """
    canvas.saveState()
    canvas.drawImage("sll_logo.gif", 2 * cm, defaultPageSize[1] - 2 * cm, 4 * cm, 1.25 * cm)
    canvas.restoreState()

def make_note(parameters):
    """Builds a pdf note based on the passed dictionary of parameters.
    For the structure of the parameters, see the code of the function
    make_example_note.
    """
    story = []
    story.append(Paragraph("Project status note", h1))
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
        if headline == 'Samples': 
            data = parameters['sample_table']
            t=Table(data,5*[1.25*inch], len(data)*[0.25*inch])
            t.setStyle(TableStyle([('ALIGN',(1,1),(-2,-2),'RIGHT'),
                       ('VALIGN',(0,0),(0,-1),'TOP'),
                       ('ALIGN',(0,-1),(-1,-1),'CENTER'),
                       ('VALIGN',(0,-1),(-1,-1),'MIDDLE'),
                       ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                       ('BOX', (0,0), (-1,-1), 0.25, colors.black),
                       ]))
            story.append(t)
        
    doc = SimpleDocTemplate("{project_name}_status_note.pdf".format(**parameters))
    doc.build(story, onFirstPage=formatted_page, onLaterPages=formatted_page)

def custom_sort(a, b):
    try:
        return cmp ( int(a), int(b) )
    except:
        try:
            int(a)
        except:
            return 1
        try: 
            int(b)
        except:
            return -1
        return cmp (a, b)
    

def make_example_note():
    """Make a note with some simple nonsensical data. Looking at this function
    and running it to make a PDF should give an idea about the structure of the
    script.
    """
    parameters = {
    "project_name": "A_test",
    "customer_reference": "Some_test",
    "uppnex_project_id": "b2013444",
    }

    make_note(parameters)

def make_status_note(prj="", opts=None):
    # Temp views for FlowcellQCMetrics and ProjectSummary while we wait for permanent ones
    fc_map_fun = '''function(doc) { if (doc.Entity_type == "FlowcellQCMetrics") emit(doc, null);}'''                       
    up_map_fun = '''function(doc) { if (doc.Entity_type == "ProjectSummary") emit(doc, null);}'''                       
    sa_map_fun = '''function(doc) { if (doc.Entity_type == "SampleQCMetrics") emit(doc, null);}'''                       
    parameters = {}
    phix_err_cutoff = 2.0
    couch = couchdb.Server("http://maggie.scilifelab.se:5984")
    qc = couch['qc']

    doc = None
    res = qc.query(up_map_fun)
    for r in res:
        if r['key']['Project_id'] == prj:
            print "DEBUG: ProjectSummary doc: (you can paste this ID into the CouchDB web interface)", r['id']
            doc = r['key']
    if not doc: 
        print "Could not find project summary for " + prj
        print "Available projects: "
        for r in res:
            print r['key']['Project_id']
        sys.exit(0)

    # print doc 

    all_passed = True
    try:
        # Project ID
        parameters['project_name'] = doc['Project_id']
        # Customer reference
        if (opts.customer_ref): parameters['customer_reference'] = opts.customer_ref
        else:
            customer_ref = "no customer reference given"
            if obj['customer_prj']: customer_ref = obj['customer_prj']
            parameters['customer_reference'] = customer_ref
        # Uppnex ID (can be manually provided by user)
        if (opts.uppnex_id): parameters['uppnex_project_id'] = opts.uppnex_id # User provided Uppnex ID
        else: parameters['uppnex_project_id'] = doc['Uppnex_id']
        # Should we list all the flow cells for this project?
        
        # A table of samples, # ordered reads, # sequenced, OK / not OK
        sample_table = []
        sample_table.append(['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
        slist = doc['Samples']
        
        for s in sorted(slist.keys(), cmp=custom_sort):
            row = []
            if slist[s].has_key('scilife_name'): row.append( slist[s]['scilife_name'])
            else: row.append("N/A")
            if slist[s].has_key('customer_name'): row.append(slist[s]['customer_name'])
            else: row.append("N/A")
            if slist[s].has_key('SampleQCMetrics'):
                # Fetch SampleQCMetrics doc(s) corresponding to sample
                ind_seqs = set()
                sq_metrics = slist[s]['SampleQCMetrics'] 
                for sdoc in sq_metrics:
                    print "DEBUG: SampleQCMetrics doc: ", sdoc
                    res_ = qc[sdoc]
                    ind_seqs.add(res_['sequence'])
                row.append(','.join(list(ind_seqs)))
            else:
                row.append("N/A")
            if doc.has_key('min_M_reads_per_sample_ordered'): 
                    try: 
                        row.append( round(float(doc['min_M_reads_per_sample_ordered']),2))
                    except:
                        row.append("N/A")
            else: row.append("N/A")
            if slist[s].has_key('M_reads_sequenced'): row.append(slist[s]['M_reads_sequenced'])
            else: row.append("N/A")
            status = 'N/A'
            # Check for status of sample.
            # Test 1: Check if it has a status attribute
            # Test 2: If not, check if it has an M_reads_sequenced attribute and a min_M_reads_per_sample_ordered attribute and compare those
            # Test 3: If not, check if there is a SampleQCMetrics attribute and go in and sum up all the reads (not implemented yet)
            if slist[s].has_key('status'): 
                status = slist[s]['status']
            elif slist[s].has_key('M_reads_sequenced') and doc.has_key('min_M_reads_per_sample_ordered'):
                    if slist[s]['M_reads_sequenced'] >= doc['min_M_reads_per_sample_ordered']: status = "P"
                    else: status = "NP"
            else:
                pass
                #print "Couldn't determine status of sample ", s
            if status != "P": all_passed = False
            row.append(status)
            sample_table.append(row)
    except: 
        print "Failed to retrieve all information. Currently looking at this data structure:", doc
        sys.exit(0)
    # How to define whether a project is finished?
    if not all_passed: 
        parameters['finished'] = 'Not finished, or cannot yet assess if finished.'
    else:
        parameters['finished'] = 'Project finished.'
    parameters['sample_table']=sample_table
    make_note(parameters)

def main():
    usage = """ Generate status note for a specific project, based on information in StatusDB. 

Usage:

python project_status_note.py <project ID (e g J.Lindberg_12_01>  

The first option is mandatory. There are further flags you can use:

-u, --uppnex: Manually insert Uppnex project ID into the report.
-o, --ordered-million-reads: Manually insert the ordered number of read pairs (in millions)
-r, --customer-reference: Manually insert customer reference (the customer's name for the project) into reports
"""

    if len(sys.argv) < 2:
        print usage
        sys.exit(0)

    parser = optparse.OptionParser()
    parser.add_option('-u', '--uppnex', action="store", dest="uppnex_id", default="", help="Manually insert UPPNEX ID into reports")
    parser.add_option('-o', '--ordered-million-reads', action="store", dest="ordered_million", default="N.A.", help="Manually insert the ordered number of read pairs (in millions) into reports")
    parser.add_option('-r', '--customer-reference', action="store", dest="customer_ref", default="N.A.", help="Manually insert customer reference (the customer's name for the project) into reports")

    (opts, args) = parser.parse_args()

    if len(sys.argv)==1:
        print "No input arguments given, will generate default example note."
        print usage
        make_example_note()
    else:
        proj = sys.argv[1]

        print "Generating status for project ", proj
        p = make_status_note(proj,opts)

if __name__ == "__main__":
    main()


