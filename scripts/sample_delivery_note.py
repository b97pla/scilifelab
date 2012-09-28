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

    doc = SimpleDocTemplate("{scilifelab_name}_note_{FC_id}.pdf".format(**parameters))
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
    
    parameters = {}
    phix_err_cutoff = 2.0
    couch = couchdb.Server("http://maggie.scilifelab.se:5984")
    qc = couch['qc']

    found_fc = False
    found_proj_for_fc = False
    fc_doc = None
    prj_doc = None

    # Identify the correct flow cell.
    for res in qc.view("_design/entitytypes/_view/FlowcellQCMetrics"):
        doc = qc[res.id]
        full_fc_name = doc['name']
        short_fc_name = full_fc_name.split("_")[-1]
        if short_fc_name == fc: 
            found_fc = True
            fc_doc = doc

    # Identify the correct ProjectSummary document. We will need this for a lot of things.
    for res in qc.view("_design/entitytypes/_view/ProjectSummary"):
        doc = qc[res.id]
        if doc['Project_id'] == prj:
            found_proj_for_fc = True
            prj_doc = doc

    avail_proj = set()
    error_rates = {}
    qvs = {}

    try:
                    
        # Now we fill in the parameters dictionary for generating the report

        # Project ID
        assert(prj_doc['Project_id'] == prj)
        parameters['project_name'] = prj_doc['Project_id']

        # Customer reference (can be manually provided by user; if not, its value is "N/A")
        if (opts.customer_ref != "N/A"): 
            parameters['customer_reference'] = opts.customer_ref
            print "DEBUG: Customer reference taken from command-line"
        else:
            print "DEBUG: Looking for customer reference in ProjectSummary document"
            try:
                customer_ref = prj_doc['Customer_reference']
            except:
                print "WARNING: Could not obtain customer reference from ProjectSummary ", prj_doc['_id']
            if len(customer_ref) > 1: 
                parameters['customer_reference'] = customer_ref
            else: 
                parameters['customer_reference'] = "No customer reference given"
                print "WARNING: Possible malformed/empty customer reference in ProjectSummary ", prj_doc['_id']

        # Uppnex ID (can be manually provided by user)
        if (opts.uppnex_id != "N/A"): 
            parameters['uppnex_project_id'] = opts.uppnex_id # User provided Uppnex ID
        else:
            print "DEBUG: Looking for UPPNEX reference in ProjectSummary document"
            try:
                customer_ref = prj_doc['Uppnex_id']
            except:
                print "WARNING: Could not obtain UPPNEX ID from ProjectSummary ", prj_doc['_id']
            if len(customer_ref) > 1: 
                parameters['uppnex_project_id'] = customer_ref
            else: 
                parameters['uppnex_project_id'] = "N/A"
                print "WARNING: Possible malformed/empty UPPNEX project ID in ProjectSummary ", prj_doc['_id']

        # Ordered amount (can be manually provided by user)
        if (opts.ordered_million != "N/A"): 
            parameters['ordered_amount'] = opts.ordered_million
        else:
            try: 
                temp = prj_doc['Min_M_reads_per_sample_ordered']
                ordered_million = str(temp)
            except:
                print "WARNING: Could not obtain M ordered reads from ProjectSummary ", prj_doc['_id']
            if len(ordered_million) > 0 and ordered_million != ' ':
                ordered_amnt = str(round(float(ordered_million),1))
            else: ordered_amnt = "N/A"    
            parameters['ordered_amount'] = ordered_amnt

        # Flowcell ID
        parameters['FC_id'] = fc # Taken from the command line.

        # Now descend into sample list from ProjectSummary and try to figure out which to include!
        slist = prj_doc['Samples']
        sufficient_yield = True
        for s in sorted(slist.keys()):
            if slist[s].has_key('SampleQCMetrics'): 
                # Check if there are SampleQCMetrics for the appropriate flow cell
                samp_metrics_entries = slist[s]['SampleQCMetrics']
                for entry in samp_metrics_entries:
                    if qc[entry]['flowcell'] == fc:
                        # Now we have the sample, flowcell and project ... start filling in specific stuff
                        parameters['success'] = "Successful run. "
                        # SciLife name
                        if slist[s].has_key('scilife_name'): parameters['scilifelab_name'] = slist[s]['scilife_name']
                        else: 
                            print "WARNING: Could not get SciLife name for sample ", s
                            parameters['scilifelab_name'] = "N/A"
                        # Customer name
                        if slist[s].has_key('customer_name'): parameters['customer_name'] = slist[s]['customer_name']
                        else: 
                            print "WARNING: Could not get customer name for sample ", s
                            parameters['customer_name'] = "N/A"
                        # Rounded read count
                        if qc[entry].has_key('bc_count'): 
                            raw = qc[entry]['bc_count']
                            rounded = str(round(float(raw)/1000000,1))
                            parameters['rounded_read_count'] = rounded
                        else:
                            print "WARNING: Could not get rounded read count for sample"
                            parameters['rounded_read_count'] = 'N/A'
                        # Sufficient yield?
                        try: 
                            dummy = float(parameters['ordered_amount'])
                            if float(parameters['rounded_read_count']) < float(parameters['ordered_amount']): sufficient_yield = False
                        except:
                             sufficient_yield = False
                        if not sufficient_yield: parameters['success'] = 'The yield may not be sufficient. '
                        # Start date
                        if qc[entry].has_key('date'): parameters['start_date'] = qc[entry]['date']
                        # PhiX error rate
                        if qc[entry].has_key('lane'):
                            lane = qc[entry]['lane']
                            try:
                                lane_metrics = fc_doc['metrics']['illumina']['Summary']
                                phix_err = 0.5 * ( float(lane_metrics['read1'][lane]['ErrRatePhiX']) + float(float(lane_metrics['read2'][lane]['ErrRatePhiX'])))
                                parameters['phix_error_rate'] = str(phix_err)
                                if phix_err > phix_err_cutoff:
                                    print "WARNING: High error rate. "
                                    parameters['success'] += 'Average PhiX error rate above 2.0%.'
                            except:
                                print "WARNING: Could not obtain lane metrics from FlowcellQCMetrics document ", fc_doc
                        else:
                            print "WARNING: Could not find lane for sample ", s
                        error_rates[s] =str( phix_err)
 
                        # Average quality score
                        try:
                            qv_counts = qc[entry]['metrics']['fastqc']['stats']['Per sequence quality scores']['Count']
                            parameters['avg_quality_score'] = calc_avg_qv(qv_counts)
                        except:
                            print "WARNING: Could not get FastQC quality scores from ", entry
                            parameters['avg_quality_score'] = "N/A"
                        qvs[s] = str(parameters['avg_quality_score'])

#                        print parameters
                        make_note(parameters)
    except:
        print "Could not fetch all info from StatusDB"
        sys.exit(0)

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
    parser.add_option('-u', '--uppnex', action="store", dest="uppnex_id", default="N/A", help="Manually insert UPPNEX ID into reports")
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


