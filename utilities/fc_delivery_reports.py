#!/usr/bin/env python
"""Make delivery notes for a flowcell

Usage:
     fc_delivery_reports.py <flowcell id> 
                       [--archive_dir=<archive directory>
                        --analysis_dir=<analysis directory>]

Given a flowcell id, make delivery reports for all projects on that flowcell.
This script relies on Illumina data being delivered to an archive directory, 
and that the bcbb pipeline has been run in the analysis directory.

The script loads the run_info.yaml and generates a report for each project.

Options:

  -a, --archive_dir=<archive dir>          The illumina data archive root in which 
                                           flowcells are found
  -b, --analysis_dir=<analysis directory>  The directory where bcbb analyses are found
                                           in FLOWCELL_ID directories
  -n, --dry_run                            Don't do anything, just list what will happen

  --v1.5                                   Use success criteria for v1.5 flow cells
"""

import os
import sys
from optparse import OptionParser
from operator import itemgetter
import yaml
import glob
import re
from mako.template import Template
from mako.lookup import TemplateLookup

from bcbio.log import create_log_handler
# from bcbio.pipeline import log
import bcbio.templates.mako2rst as m2r
from texttable import Texttable

from bcbio.google import bc_metrics 
from bcbio.solexa.flowcell import get_flowcell_info 
import read_illumina_summary_xml as summ
from bcbio.pipeline.config_loader import load_config
from bcbio.scilifelab.google.project_metadata import ProjectMetaData

def fixProjName(pname):
    newname = pname[0].upper()
    postperiod = False
    for i in range(1, len(pname)):
        if pname[i] == ".": 
            newname += pname[i]
            postperiod = True
        elif postperiod: 
            newname += pname[i].upper()
            postperiod = False
        else:
            newname += pname[i]
            postperiod = False
    return newname

TEMPLATE="""\
Delivery report for ${project_id}
=================================

${latex_opt}

Summary
--------

${summary}

Sequence yield per sample
-------------------------

${yieldtable}

"Rerun lane" means that the samples were in a 'mixed' rerun lane containing samples from different projects. In these cases the expected number of sequences will be less than the expected number for a whole lane. 

"Unmatched" gives the number of sequences that could not be reliably assigned to a specific barcode sequence (and thus to a specific sample) in a multiplexed run. There is always a certain amount of unmatched reads. In a rerun lane, the unmatched number will not be in proportion to the number of sequences for your specific samples, so the number will seem higher than it actually is.

"High" means that the number of unmatched sequences (see above) was higher than expected. 

"Low" means that the number of sequences for the sample was lower than expected for the run. Note that if this was a re-run, the total number of sequences for the sample (including previous runs) could still be sufficient for the sample to be finished. In such cases, the Summary should contain information about whether this is the case.

Delivery
--------

NOTE: This delivery note only concerns raw sequence data delivery, that is, raw FASTQ sequence files (de-multiplexed if applicable.)
If you have ordered analysis, you will be notified of the analysis results later.

The clustering was performed on a cBot cluster generation system using
a HiSeq paired-end read cluster generation kit according to the
manufacturer's instructions. The samples were sequenced on an Illumina
HiSeq 2000 as paired-end reads to 100 bp. All lanes were spiked
with 1% phiX control library, except for lane 8, which has 2% phiX.
The sequencing runs were performed according to the
manufacturer's instructions. Base conversion was done using Illumina's OLB v1.9.

Sequences derived from the phiX control library have been filtered from the delivered files. In some cases, the sequences may contain adapter 
sequences from the library preparation. This can be due to short fragment lengths
(in small-RNA-seq experiments) when it is expected, or due to "adapter contamination"
resulting from adapter or primer dimerization, which is usually not expected. We do not remove adapters by default, because this 
procedure is error-prone and better left to the customer who wants to use the data.
For mate-pair runs, the sequences will contain mate-pair linkers, and here too we
do not remove them by default because there is no reliable standard way to do so and
is therefore better left up to the end user. 
Please contact us for more information about
how to remove PhiX, adapter contamination and mate pair linkers. 

We'd like to hear from you! Please notify us when you publish using data produced at Science For Life Laboratory (SciLifeLab) Stockholm. To acknowledge SciLifeLab Stockholm in your article, you can use a sentence like "The authors would like to acknowledge support from Science for Life Laboratory, the national infrastructure SNISS, and Uppmax for providing assistance in massively parallel sequencing and computational infrastructure."

If you have any questions about this delivery, or other issues, you are welcome to email us at genomics_support@scilifelab.se.

General information
-------------------

${infotable}

The sequence files are named after the following scheme:
lane_date_flowcell-ID_sample_barcode-index_1(2).fastq, where the 1 or 2 represents the first
(forward) and the second (reverse) read in a paired-end run. Single
end runs will have only the first read. The files only contain
sequences that have passed Illumina's chastity filter. The quality scores in the fastq files
are in the "Phred64" format, sometimes known as "Illumina 1.3+" format.

Run information
---------------

Required for successful run:

- Average error rate for read1 and read2 < 2%

Summary read 1
~~~~~~~~~~~~~~

${read1table}

Summary read 2
~~~~~~~~~~~~~~

${read2table}

QC plots
~~~~~~~~

Note: The plots below are zoomable without loss of resolution. For example, in Adobe Reader, you could select View > Zoom > Zoom To > Magnification 400% to see the plots more clearly. 

.. raw:: latex
   
   \clearpage

Quality score
^^^^^^^^^^^^^
${qcplots}


.. raw:: latex
   
   \clearpage

Percentage QV>30
^^^^^^^^^^^^^^^^

${qc30plots}


.. raw:: latex
   
   \clearpage

Error rate
^^^^^^^^^^

${errorrate}

"""

def main(flowcell_id, archive_dir, analysis_dir, config_file):
    print " ".join([flowcell_id, archive_dir, analysis_dir])
    fp = os.path.join(archive_dir, flowcell_id, "run_info.yaml")
    with open(fp) as in_handle:
        run_info = yaml.load(in_handle)
    if config_file:
        config = load_config(config_file)
    else:
        config = {}
    project_ids = dict()
    for lane in run_info:
        (l, id) = [x.strip() for x in lane['description'].split(",")]
        if project_ids.has_key(id):
            if not lane in project_ids[id]: project_ids[id].append(lane)
        else:
            project_ids[id] = [lane]
        # Check here if project is a "sub project" of the lane
        if not lane.has_key('multiplex'): continue
        for s in lane['multiplex']:
            if s.has_key('sample_prj'):
                if project_ids.has_key(s['sample_prj']):
                    if lane not in project_ids[s['sample_prj']]: project_ids[s['sample_prj']].append(lane)
                else:
                    project_ids[s['sample_prj']] = [lane]
        
    sphinx_defs = []
    for k in project_ids.keys():
        lanes = [x['lane'] for x in project_ids[k]]
        proj_file_tag = k + "_" + get_flowcell_info(flowcell_id)[1] + get_flowcell_info(flowcell_id)[0][0]
        print("INFO: saw project %s in lanes %s" %( k, ", ".join(lanes)))
        sphinx_defs.append("('%s', '%s_delivery.tex', 'Raw data delivery note', u'SciLifeLab Stockholm', 'howto'),\n"  % (proj_file_tag, proj_file_tag))
        projectfile = "%s.mako" % (proj_file_tag) 
        fp = open(projectfile, "w")
        fp.write(TEMPLATE)
        fp.close()
        mylookup = TemplateLookup(directories=['./'])
        tmpl = Template(filename=projectfile, lookup=mylookup)
        proj_conf = {
            'id' : k,
            'lanes' : project_ids[k],
            'archive_dir' : archive_dir, 
            'analysis_dir' : analysis_dir,
            'flowcell' : flowcell_id,
            'config' : config,
            }
        d = generate_report(proj_conf)
        rstfile = "%s.rst" % (proj_file_tag)
        fp = open(rstfile, "w")
        fp.write(tmpl.render(**d))
        fp.close()

    sphinxconf = os.path.join(os.getcwd(), "conf.py")
    if not os.path.exists(sphinxconf):
        print("WARNING: no sphinx configuration file conf.py found: you have to edit conf.py yourself!")
    else:
        fp = open(sphinxconf)
        lines = fp.readlines()
        fp.close()
        sdout = []
        modify_conf = False
        for sd in sphinx_defs:
            if not sd in lines:
                sdout.append(sd)
                modify_conf = True
        if modify_conf:
            i = lines.index("latex_documents = [\n")
            newconf = lines[:i+3] + sdout + lines[i+3:]
            ## Change the preamble
            i = newconf.index("#'preamble': '',\n")
            newconf = newconf[:i+1] + _latex_preamble() + newconf[i+1:]
            ## Set the logo
            i = newconf.index("#latex_logo = None\n")
            newconf = newconf[:i+1] + _latex_logo() + newconf[i+1:]
            fp = open("conf.py", "w")
            fp.write("".join(newconf))
            fp.close()


def _latex_logo():
    '''Set the logo'''
    logo = ["latex_logo = '/proj/a2010002/projects/delivery_reports/grf/scilife-sniss.jpg'"]
    return logo

def _latex_preamble():
    '''Template for preamble. Sets new header'''
    preamble = ["'preamble' : r'''",
                r'''
                \usepackage[headheight=3cm]{geometry}
                \makeatletter
                \fancypagestyle{plain}{%
                \fancyhf{}
                \fancyhead[L]{{
                \begin{tabular}{l}
                Science for Life Laboratory (SciLifeLab)\\
                \textbf{Document type}\\
                BLA\\
                \textbf{Creation date}\\
                2012-
                \end{tabular}
                }}
                \fancyhead[C]{{
                \begin{tabular}{l}
                \textbf{Document name}\\
                HiSeq Delivery note\\
                \textbf{Valid from}\\
                \\
                \end{tabular}
                }}
                \fancyhead[R]{{\begin{tabular}{ll}
                \textbf{ID number\_edition} & \\
                20140\_1 & \\
                \textbf{Issuer} & \textbf{Approver}\\
                Mikael H/MH & \\
                \end{tabular}
                }}
                }
                \fancyhf{}
                \fancyhead[L]{{
                \begin{tabular}{l}
                Science for Life Laboratory (SciLifeLab)\\
                \textbf{Document type}\\
                BLA\\
                \textbf{Creation date}\\
                2012-
                \end{tabular}
                }}
                \fancyhead[C]{{
                \begin{tabular}{l}\\
                \textbf{Document name}\\
                HiSeq Delivery note\\
                \textbf{Valid from}\\
                \\
                \end{tabular}
                }}
                \fancyhead[R]{{\begin{tabular}{ll}
                \\
                \textbf{ID number\_edition} & \\
                20140\_1 & \\
                \textbf{Issuer} & \textbf{Approver}\\
                Mikael H/MH & \\
                \end{tabular}
                }}
                \makeatother
                \pagestyle{fancy}
                ''',
                "'''"] 
    return preamble

def generate_report(proj_conf):
    
    #######
    ### Metadata fetched from the 'Genomics project list' on Google Docs
    ###
    uppnex_proj = ''
    min_reads_per_sample = ''
    try:
    	proj_data = ProjectMetaData(proj_conf['id'], proj_conf['config'])
    	uppnex_proj = proj_data.uppnex_id
        project_id = proj_data.project_id
        queue_date = proj_data.queue_date
        no_samples = proj_data.no_samples
        lanes_plates = proj_data.lanes_plates
        min_reads_per_sample = proj_data.min_reads_per_sample
        customer_reference = proj_data.customer_reference
        application = proj_data.application
        no_finished_samples = proj_data.no_finished_samples
    except:
        print("WARNING: Could not fetch meta data from Google Docs")

    d = { 
        'project_id' : proj_conf['id'],
        'latex_opt' : "",
        'summary' : "",
        'infotable' : "",
        'lanetable' : "",
        'read1table': "",
        'read2table': "",
        'qcplots': "",
        'qc30plots': "",
        'errorrate': "",
        'yieldtable': "",
        }

    ## Latex option (no of floats per page)
    floats_per_page = '.. raw:: latex\n\n   \setcounter{totalnumber}{8}'
    d.update(latex_opt = floats_per_page)

    ## General info table
    tab = Texttable()
    if not uppnex_proj or len(uppnex_proj) < 4 or uppnex_proj[0:4] != 'b201':
        uppnex_proj = "b201YXXX"
        print "WARNING: Could not find UPPNEX project"

    run_name_comp = proj_conf['flowcell'].split('_')
    simple_run_name = run_name_comp[0] + "_" + run_name_comp[3]
    proj_level_dir = fixProjName(proj_conf['id'])
    instr_id = run_name_comp[1]
    fc_name, fc_date = get_flowcell_info(proj_conf['flowcell'])
    tab.add_row(["Run name:", proj_conf['flowcell']])
    del_base = "/proj/"
    proj_id = proj_conf['id']
    try: 
        if len(customer_reference) > 1:
            proj_id += ' (' + customer_reference + ')'
    except:
        pass

    tab.add_rows([["Project id:", proj_id], 
                  ["Date:", fc_date],
                  ["Instrument ID:", instr_id],
                  ["Flow cell ID:", fc_name],
                  ["Uppnex project:", uppnex_proj],
                  ["Delivery directory:", del_base + uppnex_proj + "/INBOX/" + proj_level_dir + "/" + simple_run_name]])
    d.update(infotable=tab.draw())
    
    ## Lane table
    tab = Texttable()
    tab.add_row(["Lane", "Sample(s)"])
    for l in proj_conf['lanes']:
        main_proj = l['description'].split(',')[1].strip()
        samples = []
        if l.has_key('multiplex'):
            for mp in l['multiplex']:
                if mp.has_key('sample_prj'):
                    if mp['sample_prj'] == proj_conf['id']:
                        samples.append(mp['name'])
            tab.add_row([l['lane'], ", ".join(samples)])
        else:
            tab.add_row([l['lane'], "Non-multiplexed lane"])
    d.update(lanetable=tab.draw())
    
    tab_r1 = Texttable()
    tab_r2 = Texttable()
    tab_r1.set_cols_width([2,12,12,12,12,12,12,30])
    tab_r2.set_cols_width([2,12,12,12,12,12,12,30])
    tab_r1.add_row(["Lane", "Clu. dens. #/mm2","% PF clusters","Clu. PF #/mm2", "% phas/prephas", "% aln PhiX", "% error rate", "Comment"])
    tab_r2.add_row(["Lane", "Clu. dens. #/mm2","% PF clusters","Clu. PF #/mm2", "% phas/prephas", "% aln PhiX", "% error rate", "Comment"])

    # These should be moved to a cfg file. ( + perhaps provide an alternative for v1.5 FC )
    if (options.v1_5_fc): min_clupf = 300 
    else: min_clupf = 475
    max_phas = 0.4
    max_prephas = 1.0 # 0.5
    max_mean_err = 2

    statspath = os.path.join(proj_conf['archive_dir'], proj_conf['flowcell'], "Data", "reports", "Summary")
    stats = summ.getQCstats(statspath)

    # Check quality criteria and add comments
    comm_r1 = ''
    comm_r2 = ''
    ok_r1 = True
    ok_r2 = True
    ok_cludens_r1 = True
    ok_cludens_r2 = True
    ok_err_rate = True 
    ok_err_r1 = True
    ok_err_r2 = True

    for l in proj_conf['lanes']:

        # Cluster densities
        clu_dens_r1 =  stats['raw_cluster_dens']['read1'][l['lane']]
        clu_dens_r2 =  stats['raw_cluster_dens']['read2'][l['lane']]
        clu_dens_sd_r1 =  stats['raw_cluster_dens_sd']['read1'][l['lane']]
        clu_dens_sd_r2 =  stats['raw_cluster_dens_sd']['read2'][l['lane']]
        clu_dens_string_r1 = str(clu_dens_r1) + '+/-' + str(clu_dens_sd_r1) 
        clu_dens_string_r2 = str(clu_dens_r2) + '+/-' + str(clu_dens_sd_r2) 

        # Cluster PF densities
        clu_dens_pf_r1 =  stats['pf_cluster_dens']['read1'][l['lane']]
        clu_dens_pf_r2 =  stats['pf_cluster_dens']['read2'][l['lane']]
        clu_dens_pf_sd_r1 =  stats['pf_cluster_dens_sd']['read1'][l['lane']]
        clu_dens_pf_sd_r2 =  stats['pf_cluster_dens_sd']['read2'][l['lane']]
        clu_dens_pf_string_r1 = str(clu_dens_pf_r1) + '+/-' + str(clu_dens_pf_sd_r1)
        clu_dens_pf_string_r2 = str(clu_dens_pf_r2) + '+/-' + str(clu_dens_pf_sd_r2)

        # % PF clusters
        prc_pf_r1 =  stats['prc_pf']['read1'][l['lane']]
        prc_pf_r2 =  stats['prc_pf']['read2'][l['lane']]
        prc_pf_sd_r1 =  stats['prc_pf_sd']['read1'][l['lane']]
        prc_pf_sd_r2 =  stats['prc_pf_sd']['read2'][l['lane']]
        prc_pf_string_r1 = str(prc_pf_r1) + '+/-' + str(prc_pf_sd_r1)
        prc_pf_string_r2 = str(prc_pf_r2) + '+/-' + str(prc_pf_sd_r2)

        # % phasing and prephasing
        phas_r1 = stats['phasing']['read1'][l['lane']]
        phas_r2 = stats['phasing']['read2'][l['lane']]
        prephas_r1 = stats['prephasing']['read1'][l['lane']]
        prephas_r2 = stats['prephasing']['read2'][l['lane']]
        phas_string_r1 = str(phas_r1) + '/' + str(prephas_r1)
        phas_string_r2 = str(phas_r2) + '/' + str(prephas_r2)

        # % aligned
        aln_r1 = stats['prc_aligned']['read1'][l['lane']]
        aln_r2 = stats['prc_aligned']['read2'][l['lane']]
        aln_sd_r1 = stats['prc_aligned_sd']['read1'][l['lane']]
        aln_sd_r2 = stats['prc_aligned_sd']['read2'][l['lane']]
        aln_string_r1 = str(aln_r1) + '+/-' + str(aln_sd_r1)
        aln_string_r2 = str(aln_r2) + '+/-' + str(aln_sd_r2)

        # error rate
        err_r1 = stats['error_rate']['read1'][l['lane']]
        err_r2 = stats['error_rate']['read2'][l['lane']]
        err_sd_r1 = stats['error_rate_sd']['read1'][l['lane']]
        err_sd_r2 = stats['error_rate_sd']['read2'][l['lane']]
        err_str_r1 = str(err_r1) + '+/-' + str(err_sd_r1)
        err_str_r2 = str(err_r2) + '+/-' + str(err_sd_r2)
        
        comm_r1 = ""
        comm_r2 = ""

        # check criteria
        if float(clu_dens_pf_r1[:-1]) < min_clupf: 
            ok_r1 = False
            ok_cludens_r1 = False
            comm_r1 += "Low cluster density. "
        if float(clu_dens_pf_r2[:-1]) < min_clupf: 
            ok_r2 = False
            ok_cludens_r2 = False
            comm_r2 += "Low cluster density. "
        avg_error_rate = (float(err_r1) + float(err_r2))/2
        if avg_error_rate > max_mean_err:
            ok_err_rate = False
        if float(err_r1) > max_mean_err:
            comm_r1 += "High error rate. "
            ok_err_r1 = False
        if float(err_r2) > max_mean_err:
            comm_r2 += "High error rate. "
            ok_err_r2 = False

        if comm_r1 == "": comm_r1 = "OK"        
        if comm_r2 == "": comm_r2 = "OK"

        tab_r1.add_row([l['lane'], clu_dens_string_r1, prc_pf_string_r1, clu_dens_pf_string_r1, phas_string_r1, aln_string_r1, err_str_r1, comm_r1])
        tab_r2.add_row([l['lane'], clu_dens_string_r2, prc_pf_string_r2, clu_dens_pf_string_r2, phas_string_r2, aln_string_r2, err_str_r2, comm_r2])

    # Reinitialize comments for the summary. (Which will be for several lanes, potentially)
    comm_r1 = ""
    comm_r2 = ""
 
    if not ok_cludens_r1: comm_r1 += "Low cluster density. " 
    if not ok_cludens_r2: comm_r2 += "Low cluster density. " 
    if not ok_err_rate:
        if not ok_err_r1: 
            ok_r1 = False
            comm_r1 += "High error rate. "
        if not ok_err_r2: 
            ok_r2 = False
            comm_r2 += "High error rate. "

    if (ok_r1 and ok_r2): 
        comm_r1 = comm_r2 = "OK"
        d.update(summary = "Successful run in terms of error rate. ")
    else:  
        if (ok_r1): 
            comm_r1 = "OK"
            d.update (summary = "Read 2 did not pass quality criteria: " + comm_r2)
        elif (ok_r2):
            comm_r2 = "OK"
            d.update (summary = "Read 1 did not pass quality criteria: " + comm_r1)
        else:
            d.update (summary = "Did not pass quality criteria. Read 1: " + comm_r1 + " Read 2: " + comm_r2)


    d.update(read1table=tab_r1.draw())
    d.update(read2table=tab_r2.draw())
        
    ## qcplots
    byCycleDir = os.path.join(proj_conf['archive_dir'], proj_conf['flowcell'], "Data", "reports", "ByCycle")
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "QScore_L%s.png" % (l['lane']))), width="100%"))
    d.update(qcplots= "\n".join(res))

    ## qc30plots
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "NumGT30_L%s.png" % (l['lane']))), width="100%"))
    d.update(qc30plots= "\n".join(res))

    ## qcplots
    res = []
    for l in proj_conf['lanes']:
        res.append(m2r.image(os.path.relpath(os.path.join(byCycleDir, "ErrRate_L%s.png" % (l['lane']))), width="100%"))
    d.update(errorrate= "\n".join(res))

    ## Sequence yield table
    target_yield_per_lane = 143000000.0
    if (options.v1_5_fc):  target_yield_per_lane = 60000000.0
    tab = Texttable()
    tab.add_row(['Lane','Sample','Number of sequences','Million sequences ordered','Comment'])
    
    run_info_yaml = os.path.join(proj_conf['archive_dir'],proj_conf['flowcell'],"run_info.yaml")

    if not os.path.exists(run_info_yaml):
        print("WARNING: could not find required run_info.yaml configuration file at '%s'" % run_info_yaml)
        return

    with open(run_info_yaml) as in_handle:
        run_info = yaml.load(in_handle)

    fc_name, fc_date = get_flowcell_info(proj_conf['flowcell'])
    low_yield = False
    
    bc_multiplier = 0.75 # Should move to cfg file

    ok_samples = []
    low_samples = []

    for l in proj_conf['lanes']:
        
	bc_file_name_prefix = os.path.join(proj_conf['analysis_dir'], proj_conf['flowcell'], '_'.join([l['lane'], fc_date, fc_name, "nophix_barcode"]), '_'.join([l['lane'], fc_date, fc_name, "nophix"]))
        bc_file = bc_file_name_prefix + ".bc_metrics"
	if not os.path.exists(bc_file):
		bc_file = bc_file_name_prefix + "_bc.metrics"
        try:
            bc_file = open(bc_file_name)
        except:
            sys.exit("Could not find bc metrics file " + bc_file_name)
        bc_count = {}
        for line in bc_file:
            c = line.strip().split()
            bc_count[c[0]]=c[1] + ' (~' + str (int ( round (float(c[1])/1000000) ) ) + " million)"
        no_samples = len(bc_count) - 1
        if no_samples == 0:
            print("WARNING: did not find a BC metrics file... Skipping lane %s for %s" %(l['lane'], proj_conf['id']))
            continue
        
        target_yield_per_sample = ''
        try:
            min_reads_per_sample = round(float(str(min_reads_per_sample)))
            target_yield_per_sample = min_reads_per_sample * 1000000
        except ValueError:
            min_reads_per_sample = ''
            target_yield_per_sample = bc_multiplier * target_yield_per_lane / no_samples
            
        sample_name = {}
        is_multiplexed = True
        is_rerun = False
        # Check here for each sample if it belongs to the project
        for entry in run_info:
            if entry['lane'] == l['lane']:
                projs = set()
                if entry.has_key('multiplex'):
                    for sample in entry['multiplex']:
                        if sample.has_key('sample_prj'):
                            projs.add(sample['sample_prj'])
                            if sample['sample_prj'].strip() == proj_conf['id']:
                                sample_name[sample['barcode_id']]=sample['name']
                else: is_multiplexed = False
                if len(projs) > 1: is_rerun = True
        samp_count = {}

        for k in bc_count.keys():
            if not k.isdigit(): pass
            else: 
                if sample_name.has_key(int(k)): samp_count[sample_name[int(k)]] =  bc_count[k]

        print "DEBUG: Target yield per sample = ", target_yield_per_sample
        print "DEBUG: Min reads per sample = ", min_reads_per_sample
        print "DEBUG: No samples: ", no_samples

        for k in sorted(samp_count.keys()):
            comment = ''
            if int(samp_count[k].split('(')[0]) < target_yield_per_sample: 
                comment = 'Low. '
                low_yield = True
                low_samples.append(k)
            else: ok_samples.append(k)
            if is_rerun: comment += '(rerun lane)'
            tab.add_row([l['lane'], k, samp_count[k], min_reads_per_sample, comment])
        
        if is_multiplexed:
            comment = ''
            try:
                if int (bc_count['unmatched'].split('(')[0]) > target_yield_per_sample: comment = 'High.'
                if is_rerun: comment += '(rerun lane)'
                tab.add_row([l['lane'], 'unmatched', bc_count['unmatched'], min_reads_per_sample, comment])
            except:
                print('WARNING: insufficient or no barcode metrics for lane')
        else:
            comment = ''
            for k in bc_count.keys():
                if int (bc_count[k].split('(')[0]) < bc_multiplier * target_yield_per_lane: comment = 'Low.' 
                tab.add_row([l['lane'], "Non-multiplexed lane", bc_count[k], min_reads_per_sample, comment])

    delivery_type = "Final delivery. "
    if low_yield:
        delivery_type = "Partial delivery. "
        fail_comm = "Samples " + ", ".join(low_samples) + " yielded fewer sequences than expected. These will be re-run unless this was already a re-run and the total yield is now sufficient. "
    else: fail_comm = ""

    if low_yield: 
        if len(ok_samples)>0: ok_comm = "Samples " + ", ".join(ok_samples) + " yielded the expected number of sequences or more. "
        else: ok_comm = ""
    else: ok_comm = "All samples yielded the expected number of sequences or more. "

    comm = d['summary'] + fail_comm + ok_comm
    d.update(summary = comm)

    d.update(yieldtable=tab.draw())
    return d

if __name__ == "__main__":
    usage = """
    fc_delivery_reports.py <flowcell id>
                           [--archive_dir=<archive directory> 
                            --analysis_dir=<analysis directory>]

    For more extensive help type fc_delivery_reports.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--archive_dir", dest="archive_dir", default="/bubo/proj/a2010002/archive")
    parser.add_option("-b", "--analysis_dir", dest="analysis_dir", default="/bubo/proj/a2010002/nobackup/illumina")
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",default=False)
    parser.add_option("--v1.5", dest="v1_5_fc", action="store_true", default=False)
    parser.add_option("-c", "--config-file", dest="config_file", default=None)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(
        archive_dir = os.path.normpath(options.archive_dir),
        analysis_dir = os.path.normpath(options.analysis_dir),
        config_file = options.config_file
        )
    main(*args, **kwargs)
