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
from bcbio.pipeline import log
import bcbio.templates.mako2rst as m2r
from texttable import Texttable

from bcbio.google import bc_metrics 
from bcbio.solexa.flowcell import get_flowcell_info 
import read_illumina_summary_xml as summ
from bcbio.pipeline.config_loader import load_config

TEMPLATE="""\
Delivery report for ${project_id}
=================================

Delivery
--------

The clustering was performed on a cBot cluster generation system using
a HiSeq paired-end read cluster generation kit according to the
manufacturer's instructions. The samples were sequenced on an Illumina
HiSeq 2000 as paired-end reads to 100 bp. All lanes were spiked
with 1% phiX control library, except for lane 8, which has 2% phiX.
The sequencing runs were performed according to the
manufacturer's instructions. Base conversion was done using Illumina's OLB v1.9.

Mate-pair libraries have been prepared using Roche's protocol for 454.

Note that the delivered sequences will contain sequences derived
from the PhiX control library, unless you have specifically requested
that they be removed. In some cases, the sequences may contain adapter 
sequences from the library preparation. For mate-pair runs, the sequences 
will contain mate-pair linkers. Please contact us for more information about
how to remove PhiX, adapter contamination and mate pair linkers. 

We'd like to hear from you! Please acknowledge Science for Life Laboratory (SciLifeLab Stockholm) in your articles when using data from us. Please also notify us when you publish using data produced at SciLifeLab.

Comment
--------

${comment}

General information
-------------------

${infotable}

${lanetable}

The sequence files are named after the following scheme:
lane_date_flowcell-ID_sample_barcode-index_1(2).fastq, where the 1 or 2 represents the first
(forward) and the second (reverse) read in a paired-end run. Single
end runs will have one the first read. The files only contain
sequences that have passed Illumina's chastity filter.

Run information
---------------

Required for successful run:

- Clu. PF (#/mm2) > 475K

- Phasing < 0.4%

- Prephasing < 0.5%

- Average error rate for read1 and read2 < 2%

Summary read 1
~~~~~~~~~~~~~~

${read1table}

Summary read 2
~~~~~~~~~~~~~~

${read2table}

QC plots
~~~~~~~~

Quality score
^^^^^^^^^^^^^
${qcplots}

Percentage QV>30
^^^^^^^^^^^^^^^^

${qc30plots}

Error rate
^^^^^^^^^^

${errorrate}

Sequence yield per sample
~~~~~~~~~~~~~~~~~~~~~~~~~

${yieldtable}

"""

def main(flowcell_id, archive_dir, analysis_dir):
    print " ".join([flowcell_id, archive_dir, analysis_dir])
    fp = os.path.join(archive_dir, flowcell_id, "run_info.yaml")
    with open(fp) as in_handle:
        run_info = yaml.load(in_handle)
    project_ids = dict()
    for lane in run_info:
        (l, id) = [x.strip() for x in lane['description'].split(",")]
        if project_ids.has_key(id):
            project_ids[id].append(lane)
        else:
            project_ids[id] = [lane]
        # Check here if project is a "sub project" of the lane
        for s in lane['multiplex']:
            if s.has_key('description'):
                if project_ids.has_key(s['description']):
                    project_ids[s['description']].append(lane)
                else:
                    project_ids[s['description']] = [lane]
                                                                                             
    sphinx_defs = []
    for k in project_ids.keys():
        lanes = [x['lane'] for x in project_ids[k]]
        log.info("saw project %s in lanes %s" %( k, ", ".join(lanes)))
        sphinx_defs.append("('%s', '%s_delivery.tex', 'Delivery note', u'SciLifeLab Stockholm', 'howto'),\n"  % (k, k))
        projectfile = "%s.mako" % (k)
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
            }
        d = generate_report(proj_conf)
        rstfile = "%s.rst" % (k)
        fp = open(rstfile, "w")
        fp.write(tmpl.render(**d))
        fp.close()

    sphinxconf = os.path.join(os.getcwd(), "conf.py")
    if not os.path.exists(sphinxconf):
        log.warn("no sphinx configuration file conf.py found: you have to edit conf.py yourself!")
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
            fp = open("conf.py", "w")
            fp.write("".join(newconf))
            fp.close()


def generate_report(proj_conf):
    d = { 
        'project_id' : proj_conf['id'],
        'comment' : "",
        'infotable' : "",
        'lanetable' : "",
        'read1table': "",
        'read2table': "",
        'qcplots': "",
        'qc30plots': "",
        'errorrate': "",
        'yieldtable': "",
        }

    ## General info table
    tab = Texttable()
    uppnex_proj = "b2011XXX"
    run_name_comp = proj_conf['flowcell'].split('_')
    simple_run_name = run_name_comp[0] + run_name_comp[3][0]
    tab.add_row(["Project id", proj_conf['id']])
    tab.add_rows([["Run name:", proj_conf['flowcell']],
                  ["Uppnex project:", uppnex_proj],
                  ["Delivery directory:", "/bubo/proj/" + uppnex_proj + "/INBOX/20" + simple_run_name + "_hiseq2000"]])
    d.update(infotable=tab.draw())
    
    ## Lane table
    tab = Texttable()
    tab.add_row(["Lane", "Sample(s)"])
    for l in proj_conf['lanes']:
        samples = []
        if l.has_key('multiplex'):
            for mp in l['multiplex']:
                samples.append(mp['name'])
            tab.add_row([l['lane'], ", ".join(samples)])
        else:
            tab.add_row([l['lane'], "Non-multiplexed lane"])
    d.update(lanetable=tab.draw())
    
    tab_r1 = Texttable()
    tab_r2 = Texttable()
    tab_r1.add_row(["Lane", "Clu. dens. #/mm2","% PF clusters","Clu. PF #/mm2", "% phas/prephas", "% aln PhiX", "% error rate", "Comment"])
    tab_r2.add_row(["Lane", "Clu. dens. #/mm2","% PF clusters","Clu. PF #/mm2", "% phas/prephas", "% aln PhiX", "% error rate", "Comment"])

    # These should be moved to a cfg file. ( + perhaps provide an alternative for v1.5 FC )
    min_clupf = 475
    max_phas = 0.4
    max_prephas = 0.5
    max_mean_err = 2

    statspath = os.path.join(proj_conf['archive_dir'], proj_conf['flowcell'], "Data", "reports", "Summary")
    stats = summ.getQCstats(statspath)

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
        
        # check criteria
        comm_r1 = ''
        comm_r2 = ''
        ok_r1 = True
        ok_r2 = True
        if float(clu_dens_pf_r1[:-1]) < min_clupf: 
            ok_r1 = False
            comm_r1 += "Low cluster density. "
        if float(clu_dens_pf_r2[:-1]) < min_clupf: 
            ok_r2 = False
            comm_r2 += "Low cluster density. "
        if float(phas_r1) > max_phas: 
            ok_r1 = False
            comm_r1 += "High phasing. "
        if float(phas_r2) > max_phas: 
            ok_r2 = False
            comm_r2 += "High phasing. "
        if float(prephas_r1) > max_prephas: 
            ok_r1 = False
            comm_r1 += "High prephasing. "
        if float(prephas_r2) > max_prephas: 
            ok_r2 = False
            comm_r2 += "High prephasing. "

        avg_error_rate = (float(err_r1) + float(err_r2))/2
        
        if avg_error_rate > max_mean_err:
            ok_r1 = False
            ok_r2 = False
            comm_r1 += "High error rate. "
            comm_r2 += "High error rate. "

        if (ok_r1 and ok_r2): 
            comm_r1 = comm_r2 = "OK"
            d.update(comment = "Successful run according to QC criteria.")
        else:  
            if (ok_r1): 
                comm_r1 = "OK"
                d.update (comment = "Read 2 did not pass quality criteria: " + comm_r2)
            elif (ok_r2):
                comm_r2 = "OK"
                d.update (comment = "Read 1 did not pass quality criteria: " + comm_r1)
            else:
                d.update (comment = "Did not pass quality criteria. Read 1: " + comm_r1 + " Read 2: " + comm_r2)

        tab_r1.add_row([l['lane'], clu_dens_string_r1, prc_pf_string_r1, clu_dens_pf_string_r1, phas_string_r1, aln_string_r1, err_str_r1, comm_r1])
        tab_r2.add_row([l['lane'], clu_dens_string_r2, prc_pf_string_r2, clu_dens_pf_string_r2, phas_string_r2, aln_string_r2, err_str_r2, comm_r2])

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
    tab = Texttable()
    tab.add_row(['Lane','Sample','Number of sequences','Comment'])
    
    run_info_yaml = os.path.join(proj_conf['archive_dir'],proj_conf['flowcell'],"run_info.yaml")

    if not os.path.exists(run_info_yaml):
        log.warn("Could not find required run_info.yaml configuration file at '%s'" % run_info_yaml)
        return

    #with open(run_info_yaml) as in_handle:
    #    run_info = {'details': yaml.load(in_handle)}

    with open(run_info_yaml) as in_handle:
        run_info = yaml.load(in_handle)

    # fc_name, fc_date = get_flowcell_info(proj_conf['flowcell'])
    # bc_yield = bc_metrics.get_bc_stats(fc_date,fc_name,proj_conf['analysis_dir'], run_info)
   
    fc_name, fc_date = get_flowcell_info(proj_conf['flowcell'])
    low_yield = False
    
    bc_multiplier = 0.75 # Should move to cfg file

    for l in proj_conf['lanes']:
        bc_file_name = os.path.join(proj_conf['analysis_dir'], proj_conf['flowcell'], '_'.join([l['lane'], fc_date, fc_name, "barcode"]), '_'.join([l['lane'], fc_date, fc_name, "bc.metrics"]))
        try:
            bc_file = open(bc_file_name)
        except:
            sys.exit("Could not find bc metrics file " + bc_file_name)
        bc_count = {}
        for line in bc_file:
            c = line.strip().split()
            bc_count[c[0]]=c[1] + ' (~' + str (int ( round (float(c[1])/1000000) ) ) + " million)"
        no_samples = len(bc_count)
        target_yield_per_sample = bc_multiplier * target_yield_per_lane / no_samples
        sample_name = {}
        is_multiplexed = True
        # Check here for each sample if it belongs to the project
        for entry in run_info:
            if entry['lane'] == l['lane']:
                is_main_proj = False       
                if entry['description'].split(',')[1].strip() == proj_conf['id']:
                    is_main_proj = True
                if entry.has_key('multiplex'):
                    for sample in entry['multiplex']:
                        if sample.has_key('description'):
                            if is_main_proj: log.info('Skipping sample ' + sample['name'] + ' in lane ' + l['lane'])
                            else:
                                if sample['description'].strip() == proj_conf['id']:
                                    sample_name[sample['barcode_id']]=sample['name']
                        elif is_main_proj: 
                            sample_name[sample['barcode_id']]=sample['name']
                else: is_multiplexed = False
        samp_count = {}

        for k in bc_count.keys():
            if not k.isdigit(): pass
            else: 
                if sample_name.has_key(int(k)): samp_count[sample_name[int(k)]] =  bc_count[k]
        for k in sorted(samp_count.keys()):
            comment = ''
            if int(samp_count[k].split('(')[0]) < target_yield_per_sample: 
                comment = 'Low'
                low_yield = True
            tab.add_row([l['lane'], k, samp_count[k], comment])
        
        if is_multiplexed:
            try:
                comment = ''
                if int (bc_count['unmatched'].split('(')[0]) > target_yield_per_sample: comment = 'High'
                tab.add_row([l['lane'], 'unmatched', bc_count['unmatched'], comment])
            except:
                log.warning('Unsufficient or no barcode metrics for lane')
        else:
            for k in bc_count.keys():
                comment = ''
                if int (bc_count[k].split('(')[0]) < bc_multiplier * target_yield_per_lane: comment = 'Low' 
                tab.add_row([l['lane'], "Non-multiplexed lane", bc_count[k], comment])

    if low_yield:
        comm = d['comment'] +  " Some samples had low yields."
        d.update(comment = comm)
    
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
    parser.add_option("-b", "--analysis_dir", dest="analysis_dir", default="/bubo/proj/a2010002/nobackup/romanvg")
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        print __doc__
        sys.exit()
    kwargs = dict(
        archive_dir = os.path.normpath(options.archive_dir),
        analysis_dir = os.path.normpath(options.analysis_dir)
        )
    main(*args, **kwargs)
