"""Script for generating a PDF delivery note for a sample.
The script uses reportlab for PDF, for which a user guide can be found here:
http://www.reportlab.com/software/opensource/rl-toolkit/guide/
"""

import sys
import couchdb
import optparse

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

# This is basically a template of ('headline', 'content') pairs.
# Use dictionary based python string formatting to fill out variables.
paragraphs = OrderedDict()
paragraphs["Project name"] = "{project_name} ({customer_reference})"

paragraphs["UPPNEX project id"] = "{uppnex_project_id}"

paragraphs["Flow cell id"] = "{FC_id}"

paragraphs["Sequence data directory"] = \
"/proj/{uppnex_project_id}/INBOX/{project_name}/{scilifelab_name}/{start_date}_{FC_id}"

paragraphs["Sample"] = "{scilifelab_name} / {customer_name}." \
                       "Ordered amount: {ordered_amount} million paired reads."

paragraphs["Method"] = "Clustered on cBot and sequenced on HiSeq 2000 " \
                       "according to manufacturer's instructions. Base " \
                       "conversion using OLB v1.9, demultiplexed and " \
                       "converted to fastq using CASAVA v1.8. The quality scale " \
                       "is Sanger / phred33 / Illumina 1.8+."

paragraphs["Results"] = "{rounded_read_count} million reads in lane with PhiX " \
                        "error rate {phix_error_rate}%. Average quality score " \
                        "{avg_quality_score}."

paragraphs["Comments"] = "{success}"

# Move Information section to Project Status Report

#paragraphs["Information"] = OrderedDict()

#paragraphs["Information"]["Acknowledgement"] = \
#"Please notify us when you publish using data produced at Science For Life " \
#"Laboratory (SciLifeLab) Stockholm. To acknowledge SciLifeLab Stockholm in your " \
#'article, you can use a sentence like "The authors would like to acknowledge ' \
#"support from Science for Life Laboratory, the national infrastructure SNISS, " \
#"and Uppmax for providing assistance in massively parallel sequencing and " \
#'computational infrastructure."'

#paragraphs["Information"]["Naming conventions"] = \
#"The data is delivered in fastq format using Illumina 1.8 quality scores. " \
#"There will be one file for the forward reads and one file for the reverse " \
#"reads. More information on our naming conventions can be found here."

#paragraphs["Information"]["Data access at UPPMAX"] = \
#"Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next " \
#"Generation sequence Cluster & Storage, www.uppmax.uu.se), from which the " \
#"user can access it. If you have problems to access your data, please contact " \
#"SciLifeLab genomics_support@scilifelab.se. If you have questions regarding " \
#"UPPNEX, please contact support@uppmax.uu.se. Information on how to access your " \
#"data can be found here."


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

def make_example_note():
    """Make a note with some simple nonsensical data. Looking at this function
and running it to make a PDF should give an idea about the structure of the
script.
"""
    parameters = {
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

    make_note(parameters)

def make_notes_for_fc_proj(fc="BC0HYUACXX", prj="J.Lindberg_12_01", opts=None):
    # Warning: ugliness ahead! :-)
    # Temp views for FlowcellQCMetrics and ProjectSummary while we wait for permanent ones
    fc_map_fun = '''function(doc) { if (doc.entity_type == "FlowcellQCMetrics") emit(doc, null);}'''
    ps_map_fun = '''function(doc) { if (doc.Entity_type == "ProjectSummary") emit(doc, null);}'''
    s_map_fun = '''function(doc) { if (doc.entity_type == "SampleQCMetrics") emit(doc, null);}'''
    
    parameters = {}
    phix_err_cutoff = 2.0
    couch = couchdb.Server("http://maggie.scilifelab.se:5984")
    qc = couch['qc']
    found_fc = False
    found_proj_for_fc = False
    results = qc.query(s_map_fun)
    avail_proj = set()
    error_rates = {}
    qvs = {}
    for r in results:
        obj = r['key']
        if obj.has_key('sample_prj'):
            if obj['flowcell'] == fc:
                found_fc = True
                avail_proj.add(obj['sample_prj'])
            if obj['sample_prj'] == prj and obj['flowcell'] == fc: # We have the correct flow cell and project ID
                found_proj_for_fc = True
                print "SampleQCMetrics document ID: ", r['id']

                #if found_fc: print "DEBUG: Found flow cell document"
                #if found_proj_for_fc: print "DEBUG: Found project"

                try:
                    
                    # Now we fill in the parameters dictionary for generating the report
                    
                    # Project ID
                    parameters['project_name'] = obj['sample_prj']

                    # Customer reference (can be manually provided by user)
                    if (opts.customer_ref != "N/A"):
                        parameters['customer_reference'] = opts.customer_ref
                        print "DEBUG: Customer reference taken from command-line"
                    else:
                        print "DEBUG: Looking for customer reference in ProjectSummary document"
                        customer_ref = "no customer reference given"
                        results = qc.query(ps_map_fun)
                        for r in results:
                            if prj == r['key']['Project_id']:
                                cust_ref = r['key']['Customer_reference']
                                parameters['customer_reference'] = cust_ref
                    print "Customer reference value: ", parameters['customer_reference']

                    # Uppnex ID (can be manually provided by user)
                    #print "Uppmax ID from command line: ", opts.uppnex_id
                    if (opts.uppnex_id): parameters['uppnex_project_id'] = opts.uppnex_id # User provided Uppnex ID
                    else:
                        uppnex = "NA"
                        results = qc.query(ps_map_fun)
                        for r in results:
                            if prj == r['key']['Project_id']:
                                uppnex = r['key']['Uppnex_id']
                                parameters['uppnex_project_id'] = uppnex

                    # Start date
                    parameters['start_date'] = obj['date']

                    # Ordered amount (can be manually provided by user)
                    if (opts.ordered_million != "N/A"):
                        parameters['ordered_amount'] = opts.ordered_million
                    else:
                        ordered_amnt = "N/A"
                        results = qc.query(ps_map_fun)
                        for r in results:
                            if prj == r['key']['Project_id']:
                                print "ProjectSummary document ID: ", r['id']
                                temp = r['key']['Min_M_reads_per_sample_ordered']
                                ordered_million = str(temp)
                                if len(ordered_million) > 0 and ordered_million != ' ':
                                    ordered_amnt = ordered_million
                        parameters['ordered_amount'] = ordered_amnt

                    # Flowcell ID
                    parameters['FC_id'] = obj['flowcell']
                    # Scilife sample name
                    parameters['scilifelab_name'] = obj['barcode_name']

                    # Customer sample name
                    try:
                        cust_name = "no customer sample name given"
                        if (obj['customer_sample_name']):
                            cust_name = obj['customer_sample_name']
                            print " DEBUG --- found customer sample name in SampleQCMetrics: ", cust_name
                        else:
                            print " DEBUG --- found no customer sample name corresponding to SciLife ID ", parameters['scilifelab_name'], " in SampleQCMetrics; trying ProjectSummary document ..."
                            results = qc.query(ps_map_fun)
                            for r in results:
                                if prj == r['key']['Project_id']:
                                    slist = r['key']['Samples']
                                    for s in slist.keys():
                                        without_index = "_".join(parameters['scilifelab_name'].split("_")[0:-1])
                                        if slist[s].has_key('scilife_name'):
                                            if slist[s]['scilife_name'] == without_index:
                                                if slist[s].has_key('customer_name'):
                                                    cust_name = slist[s]['customer_name']
                                                    print "Found customer sample name ", cust_name, " corresponding to ", slist[s]['scilife_name'], " or ", parameters['scilifelab_name']
                                                else: print "WARNING: Did not find customer sample name corresponding to ", slist[s]['scilife_name']," in ProjectSummary document"
                        parameters['customer_name'] = cust_name
                    except:
                        print "Failed to set customer sample name"
                    if cust_name == "no customer sample name given": print "WARNING: Could not set customer sample name for ", parameters['scilifelab_name']
                    # (Rounded) amount of read pairs, in millions
                    try:
                        rounded_amnt = round(float(obj['bc_count'])/1000000,1)
                        parameters['rounded_read_count'] = str(rounded_amnt)
                    except:
                        parameters['rounded_read_count'] = 'N/A'

                    print "After customer sample name"
                    # print "DEBUG: ", parameters
                    # Lane
                    lane = obj['lane']
                    # Average PhiX error rate and QV>30 (latter not implemented yet)
                    phix_avg = 'N/A'
                    qv_30_avg = 'N/A'
                    results = qc.query(fc_map_fun)
                    for r in results:
                        full_fc_name = r['key']['name']
                        short_fc_name = full_fc_name.split("_")[-1]
                        if short_fc_name == fc:
                            print "FlowCellQCMetrics ID: ", r['id']
                            #print "DEBUG: Found FlowCellQCMetrics document corresponding to SampleQCMetrics flowcell (whew)"
                            phix_r1 = float(r['key']['metrics']['illumina']['Summary']['read1'][lane]['ErrRatePhiX'])
                            phix_r2 = float(r['key']['metrics']['illumina']['Summary']['read3'][lane]['ErrRatePhiX'])
                            phix_avg = (phix_r1 + phix_r2)/2
                    parameters['phix_error_rate'] = str(phix_avg)
                    error_rates[parameters['scilifelab_name']] = parameters['phix_error_rate']
                    # Average QV
                    # print "DEBUG: ", obj['metrics']['fastqc']
                    avg_qv = "N/A"
                    try:
                        avg_qv = calc_avg_qv(obj['metrics']['fastqc']['stats']['Per sequence quality scores']['Count'])
                    except:
                        pass
                    parameters['avg_quality_score'] = str(avg_qv)
                    qvs[parameters['scilifelab_name']] = parameters['avg_quality_score']
                except:
                    sys.exit("Could not fetch all info from StatusDB")
                print "Making note for sample ", parameters['scilifelab_name'], " on flowcell ", parameters['FC_id'], " lane ", obj['lane']
                # Successful run?
                success_message = ''
                try:
                    if float(parameters['phix_error_rate']) < phix_err_cutoff and float(parameters['rounded_read_count']) > float(parameters['ordered_amount']):
                        success_message += "Successful run."
                    else:
                        if float(parameters['phix_error_rate']) > phix_err_cutoff: success_message += "High average error rate."
                        if float(parameters['rounded_read_count']) < float(parameters['ordered_amount']): success_message += "The yield may be lower than expected."
                except:
                    print "Warning: Could not assess success of run."
                    #print obj['_id']
                    success_message = "Could not assess success or failure of run."
                parameters['success'] = success_message
                make_note(parameters)

    print "*** Quality stats ***"
    print "Scilifelab ID\tPhiXError\tAvgQV"
    for k in sorted(error_rates.keys()):
        print k + "\t" + error_rates[k] + "\t" + qvs[k]

    if not found_fc:
        print "Could not find specified flow cell!"
        print "Available as FlowcellQCMetrics documents:"
        res = qc.query(fc_map_fun)
        fcqm = []
        for r in res:
            fcqm.append( r['key']['name'])
        if len(fcqm) != len(set(fcqm)): print "(Warning: possible duplicates)"
        for f in sorted(fcqm): print f
        print "Available as fields in SampleQCMetrics documents:"
        res = qc.query(s_map_fun)
        sqm = set()
        for r in res:
            sqm.add(r['key']['flowcell'])
        for f in sorted(sqm): print f
    elif not found_proj_for_fc:
        print "Could not find specified project in this flow cell. Available projects for flow cell:"
        for i in avail_proj: print i
            
def main():
    usage = """ Generate sample-centered delivery notes for a specific run, in a specific project, based on information in StatusDB.

Usage:

python sample_delivery_note.py <project ID (e g J.Lindberg_12_01> <flow cell ID (e g BC0HYUACXX)>

The two first options are mandatory (kind of ... the program will just generate an example note if they are omitted). There are further flags you can use:

-u, --uppnex: Manually insert Uppnex project ID into the report.
-o, --ordered-million-reads: Manually insert the ordered number of read pairs (in millions)
-r, --customer-reference: Manually insert customer reference (the customer's name for the project) into reports
"""

    if len(sys.argv) < 3 and len(sys.argv) > 1:
        print usage
        sys.exit(0)

    parser = optparse.OptionParser()
    parser.add_option('-u', '--uppnex', action="store", dest="uppnex_id", default="", help="Manually insert UPPNEX ID into reports")
    parser.add_option('-o', '--ordered-million-reads', action="store", dest="ordered_million", default="N/A", help="Manually insert the ordered number of read pairs (in millions) into reports")
    parser.add_option('-r', '--customer-reference', action="store", dest="customer_ref", default="N/A", help="Manually insert customer reference (the customer's name for the project) into reports")

    (opts, args) = parser.parse_args()

    if len(sys.argv)==1:
        print "No input arguments given, will generate default example note."
        print usage
        make_example_note()
    else:
        proj = sys.argv[1]
        fc = sys.argv[2]

        if len(fc) > 10:
            fc = fc[-10:]
            print "Flow cell ID not in expected format, I'll take a chance and use the 10 last characters: ", fc
        elif len(fc) < 10:
            sys.exit("Flow cell ID too short, exiting")

        print "Attempting to generate notes for each sample in project ", proj, " for flow cell ", fc
        p = make_notes_for_fc_proj(fc, proj,opts)

if __name__ == "__main__":
    main()

