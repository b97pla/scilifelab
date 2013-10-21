"""Module delivery_notes - code for generating delivery reports and notes"""
import os
import re
import itertools
import ast
import json
import math
import csv
import yaml
import operator
import texttable
import datetime
from cStringIO import StringIO
from collections import Counter
from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection, calc_avg_qv
from scilifelab.utils.misc import query_ok, query_yes_no
from scilifelab.report import sequencing_success
from scilifelab.report.rst import make_sample_rest_notes, make_rest_note, render_rest_note, write_rest_note
from scilifelab.report.rl import make_note, concatenate_notes, sample_note_paragraphs, sample_note_headers, project_note_paragraphs, project_note_headers, make_sample_table
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

# Software versions used in data production. Instrument specific?
software_versions = {
    'basecall_software': 'RTA',
    'demultiplex_software' : 'bcl2fastq',
    }

def _parse_instrument_config(cfile):
    """Parse a supplied yaml file with instrument ids and associated metadata and return a list of dicts
    """
    if cfile is None or not os.path.exists(cfile):
        LOG.warn("No instrument config file supplied, will use default value")
        return [{'instrument_id': 'default', 'instrument_alias': 'NN', 'instrument_version': 'NN'}]

    with open(cfile) as fh:
        return yaml.load(fh)

# http://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
def _round_read_count_in_millions(n):
    """Round absolute read counts to million reads"""
    LOG.debug("Rounding read count: got {}".format(n))
    if n is None:
        return None
    if n == 0:
        return 0
    round_factor = [2,2,1]
    millidx = max(0, min(len(round_factor) - 1, int(math.floor(math.log10(abs(int(n)))/3.0))))
    return round(float(n)/10**(6),round_factor[millidx])

def _get_ordered_million_reads(sample_name, ordered_million_reads):
    """Retrieve ordered million reads for sample

    :param sample_name: sample name (possibly barcode name)
    :param ordered_million_reads: parsed option passed to application

    :returns: ordered number of reads or None"""
    if isinstance(ordered_million_reads, dict):
        if sample_name in ordered_million_reads:
            return ordered_million_reads[sample_name]
        else:
            return ordered_million_reads.get("default", -1)
    else:
        return ordered_million_reads

def _get_phix_error_rate(lane, phix):
    """Set phix error rate for a sample based on lane

    :param lane: lane
    :param phix: parsed option passed to application

    :returns: phix error rate or None"""
    if isinstance(phix, dict):
        if int(lane) in phix:
            return phix[int(lane)]
        else:
            return -1
    else:
        return phix

def _get_bc_count(sample_name, bc_count, sample_run):
    """Retrieve barcode count for a sample

    :param sample_name: sample name
    :param bc_count: parsed option passed to application
    :param sample_run: sample run object

    :returns: barcode count or None"""
    if isinstance(bc_count, dict):
        if sample_name in bc_count:
            return bc_count[sample_name]
        else:
            return bc_count.get("default", sample_run.get("bc_count", -1))
    else:
        return bc_count


def _assert_flowcell_format(flowcell):
    """Assert name of flowcell: "[A-Z0-9\-]+"

    :param flowcell: flowcell id

    :returns: boolean
    """
    if flowcell is None:
        # Can this really be right?!?
        return True
    if not re.match("[A-Z0-9\-]+$", flowcell):
        return False
    return True

def _set_sample_run_list(project_name, flowcell, project_alias, s_con):
    """Set sample run list.

    :param project_name: project name
    :param flowcell: flowcell id
    :param project_alias: project alias argument passed to pm
    :param s_con: sample run connection

    :returns: sample_run_list
    """
    sample_run_list = s_con.get_samples(sample_prj=project_name, fc_id=flowcell)
    if not project_alias:
        return sample_run_list
    project_alias = ast.literal_eval(project_alias)
    for p_alias in project_alias:
        sample_run_list_tmp = s_con.get_samples(sample_prj=p_alias, fc_id=flowcell)
        if sample_run_list_tmp:
            sample_run_list.extend(sample_run_list_tmp)
    return sample_run_list

def _literal_eval_option(option, default=None):
    """Literally evaluate passed option.

    :param option: option passed to pm, which could be a file name
    :param default: default value of option

    :returns: parsed option
    """
    if not option:
        return default
    if os.path.exists(option):
        with open(option) as fh:
            option = json.load(fh)
    else:
        option = ast.literal_eval(option)
    return option

def _update_sample_output_data(output_data, cutoffs):
    """Update sample output data dictionary.

    :param output_data: output data dictionary
    :param cutoffs: cutoffs dictionary

    :returns: updated output data dictionary
    """
    output_data["stdout"].write("\nQuality stats\n")
    output_data["stdout"].write("************************\n")
    output_data["stdout"].write("PhiX error cutoff: > {:3}\n".format(cutoffs['phix_err_cutoff']))
    output_data["stdout"].write("QV cutoff        : < {:3}\n".format(cutoffs['qv_cutoff']))
    output_data["stdout"].write("************************\n\n")
    output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format("Scilifelab ID", "Lane", "PhiXError", "ErrorStatus", "AvgQV", "QVStatus"))
    output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format("=============", "====", "=========", "===========", "=====", "========"))
    return output_data

def _set_project_sample_dict(project_sample_item, source):
    """Set a project sample dict, mapping a project sample to sample run metrics if present in project summary.

    :param project_sample_item: a project sample item

    :returns: project_sample_d or empty dict
    """
    project_sample_d = {}

    #The structure of the database has changed for projects opened after July 1st
    #2013 (document 10294_01 for more details)
    if source == 'lims':
        # FIXME: This section does not return the correctly formatted results
        LOG.debug("This project has LIMS as source of information")
        if "library_prep" in project_sample_item.keys():
            sample_run_metrics = {k:v.get("sample_run_metrics", {}) for k,v in \
                                    project_sample_item["library_prep"].iteritems()}
            project_sample_d = {}
            for fc in sample_run_metrics.items():
                fc, metrics = fc
                for k, v in metrics.iteritems():
                    sample_run_metrics = v.get('sample_run_metrics_id', '')
                    if sample_run_metrics:
                        project_sample_d[k] = v['sample_run_metrics_id']
                    else:
                        LOG.warn("No sample_run_metrics information for sample '{}'".format(project_sample_item))
        else:
            sample_run_metrics = project_sample_item.get("sample_run_metrics", {})
            project_sample_d = {metrics[0]:metrics[1]['sample_run_metrics_id'] \
                                    for metrics in sample_run_metrics.items()}
            if not project_sample_item.get("sample_run_metrics", {}):
                LOG.warn("No sample_run_metrics information for sample '{}'".format(project_sample_item))
    else:
        if "library_prep" in project_sample_item.keys():
            project_sample_d = {x:y for d in [v.get("sample_run_metrics", {}) \
                    for k,v in project_sample_item["library_prep"].iteritems()] \
                        for x,y in d.iteritems()}
        else:
            project_sample_d = {x:y for x,y in project_sample_item.get("sample_run_metrics", {}).iteritems()}
            if not project_sample_item.get("sample_run_metrics", {}):
                LOG.warn("No sample_run_metrics information for sample '{}'".format(project_sample_item))
    return project_sample_d


def _collect_status_note_data(**kw):
    """Gather the data needed for status notes
    """
    
    # Connect to databases
    s_con = SampleRunMetricsConnection(dbname=kw.get("samplesdb"), **kw)
    fc_con = FlowcellRunMetricsConnection(dbname=kw.get("flowcelldb"), **kw)
    p_con = ProjectSummaryConnection(dbname=kw.get("projectdb"), **kw)

    project_name = kw.get("project_name")
    project = p_con.get_entry(project_name)
    if not project:
        LOG.warn("No such project '{}'".format(project_name))
        return []
    
    prj_dict = dict()
    # Update the dict with project information
    prj_dict.update(_get_project_info(project, p_con))
    
    # Get uppnex id, possible overridden on the command line
    if kw.get("uppnex_id"):
        prj_dict["uppnex_project_id"] = kw.get("uppnex_id")
        
    # Get customer reference, possible overridden on the command line        
    if kw.get("customer_reference"):
        prj_dict["customer_reference"] = kw.get("customer_reference")

    # Get the list of sample_runs, possibly restricted by flowcell
    flowcell = kw.get("flowcell")
    sample_run_list = _set_sample_run_list(project_name,flowcell,kw.get("project_alias"),s_con)
    if len(sample_run_list) == 0:
        LOG.warn("No samples for project '{}'{}. Maybe there are no sample run metrics in statusdb?".format(project_name,                                                                                                        ', flowcell {}'.format(flowcell) if flowcell else '')) 
    sample_dict = _get_sample_info(sample_run_list, prj_dict, fc_con, **kw)
    
    # Insert the actual sample run information objects into the project sample list
    prj_samples = _get_project_sample_info(project)
    for prj_sample in prj_samples:
        name = prj_sample.get('scilifeid')
        if not name:
            continue
        for sample in sample_dict:
            if sample.get('scilifeid') != name:
                continue
            id = sample.get('_id')
            if not id:
                break
            for run in [r for prep in prj_sample.get('library_preps',[]) for r in prep.get('sequencing_runs',[])]:
                if id != run.get('sample_run_id'):
                    continue
                run['sample_run'] = sample
    
    fc_dict = []
    if flowcell:
        # Get flowcell
        fc_doc = fc_con.get_flowcell_by_id(flowcell)
        if not fc_doc:
            LOG.warn("No such flowcell '{}'".format(flowcell))
            return []
        fc = fc_doc.get("name")
        
        # Get a dict with the flowcell level information
        param = _get_flowcell_info(fc_con,fc,project_name)
        # Update the dict with the relevant QC thresholds
        param.update(_get_qc_thresholds(param,kw.get("config")))
        fc_dict.append(param)
    else:
        # Get all flowcells and lanes from the sample run list
        fcs = {}
        for sample in sample_run_list:
            fc = '{}_{}'.format(sample.get("date"),
                                sample.get("flowcell"))
            if not fc in fcs:
                # Get a dict with the flowcell level information
                param = _get_flowcell_info(fc_con,fc,project_name)
                # Update the dict with the relevant QC thresholds
                param.update(_get_qc_thresholds(param,kw.get("config")))
                fc_dict.append(param)
                fcs[fc] = []
                
            fcs[fc] = list(set(fcs[fc] + [sample.get("lane")]))
        
    return prj_dict, fc_dict, sample_dict, prj_samples
    
def sample_summary_notes(project_name=None, config=None, **kw):
    """Create one sample_summary_note per sample and concatenate them at the end
    """
    #import ipdb; ipdb.set_trace()
    data = _collect_status_note_data(project_name=project_name, config=config, **kw)
    if len(data) == 0:
        return output_data
    
    prj_dict, fc_dict, _, samples = data[0:4]
    
    prep_header = ['Prep','Barcode','Average fragment size','Started','Finished','QC status']
    runs_header = ['Prep','Barcode','Flowcell','Lane','Started','Finished','QC status']
    
    restfile = '{}_sample_summary_{}.rst'.format(project_name,
                                                 datetime.datetime.now().strftime('%Y%m%d'))
    # Create a summary note per sample
    rest_notes = []
    for sample in samples:
        
        rst_dict = {}
        for field in ['scilifeid','customerid','incoming_qc_start_date','incoming_qc_finish_date','incoming_qc_status']:
            rst_dict[field] = sample[field]
        
        # Prepare the prep and sequencing table
        prep_body = []
        runs_body = []
        for prep in sample['library_preps']:
            lbl = prep['label'] or 'N/A'
            fsize = prep['fragment_size'] or sample['fragment_size_customer'] or 'N/A'
            barcode = prep['reagent_label'] or 'N/A'
            status = 'N/A' if prep['prep_qc_status'] is None else '{}'.format('Passed' if prep['prep_qc_status'] else 'Failed')
            dates = []
            for key in ['prep_start_date','prep_qc_finish_date']:
                d = 'N/A'
                try:
                    d = prep[key].strftime('%Y-%m-%d')
                except:
                    pass
                dates.append(d)
            prep_body.append([lbl,barcode,fsize,dates[0],dates[1],status])
            
            for run in prep['sequencing_runs']:
                srun = run.get('sample_run',{})
                index = srun.get('barcode','')
                fcid = srun.get('flowcell',' ')[1:]
                lane = srun.get('lane','N/A')
                
                started = 'N/A'
                try:
                    started = datetime.datetime.strptime(srun.get('date'),'%y%m%d').strftime('%Y-%m-%d')
                except:
                    pass
                
                finished = 'N/A'
                try:
                    finished = run['run_finish_date'].strftime('%Y-%m-%d')
                except:
                    pass
                
                phix = srun.get('phix_error_rate','N/A')
                q30 = srun.get('q30','N/A')
                avgq = srun.get('avgQ','N/A')
                yld = srun.get('read_count','N/A')
                status = run.get('sequencing_qc_status','N/A')
                runs_body.append([lbl,index,fcid,lane,started,finished,status])
                
        prep_table = [prep_header] + sorted(prep_body, key=operator.itemgetter(0))
        runs_table = [runs_header] + sorted(runs_body, key=operator.itemgetter(0,2,3))
        
        # Render the note
        rest_notes.append(render_rest_note(tables={'prep_table': prep_table, 'sequencing_table': runs_table, 'delivery_table': []}, 
                                           report="sample_report", 
                                           **rst_dict))
    write_rest_note(restfile,contents=rest_notes)
        
     
def sample_status_note(project_name=None, flowcell=None, username=None, password=None, url=None,
                       ordered_million_reads=None, uppnex_id=None, customer_reference=None, bc_count=None,
                       project_alias=[], projectdb="projects", samplesdb="samples", flowcelldb="flowcells",
                       phix=None, is_paired=True, config=None, **kw):
    """Make a sample status note. Used keywords:

    :param project_name: project name
    :param flowcell: flowcell id
    :param username: db username
    :param password: db password
    :param url: db url
    :param ordered_million_reads: number of ordered reads in millions
    :param uppnex_id: the uppnex id
    :param customer_reference: customer project name
    :param project_alias: project alias name
    :param phix: phix error rate
    :param is_paired: True if run is paired-end, False for single-end
    """
    
    kw["project_name"] = project_name
    kw["flowcell"] = flowcell
    kw["username"] = username
    kw["password"] = password
    kw["url"] = url
    kw["uppnex_id"] = uppnex_id
    kw["customer_reference"] = customer_reference
    kw["project_alias"] = project_alias
    kw["projectdb"] = projectdb
    kw["samplesdb"] = samplesdb
    kw["flowcelldb"] = flowcelldb
    kw["phix"] = phix
    kw["config"] = config
    
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    if not _assert_flowcell_format(flowcell):
        LOG.warn("Wrong flowcell format {}; skipping. Please use the flowcell id (format \"[A-Z0-9\-]+\")".format(flowcell) )
        return output_data
    #output_data = _update_sample_output_data(output_data, cutoffs)

    data = _collect_status_note_data(**kw)
    if len(data) == 0:
        return output_data
    
    prj_dict, fc_dict, sample_dict = data[0:3]
    
    # FIXME: parse this from configs
    for fc_param in fc_dict:
        fc_param.update(software_versions)
        fc_param["basecaller_version"] = fc_param["software"].get(fc_param["basecall_software"])
        fc_param["demultiplex_version"] = fc_param["software"].get(fc_param["demultiplex_software"],"1.8.3")
    
    # Verify that the same sample does not have two entries for the same flowcell, lane and barcode. If so, prompt for input.
    sample_d = {}
    for sample in sample_dict:
        key = "{}_{}_{}_{}".format(sample.get("date"),
                                   sample.get("flowcell"),
                                   sample.get("lane"),
                                   sample.get("barcode"))
        if key not in sample_d:
            sample_d[key] = []
        
        sample_d[key].append(sample)
        
    for key, samples in sample_d.items():
        # If we have just one sample for this key, everything is in order
        if len(samples) == 1:
            continue
            
        # Else, we need to resolve the conflict
        LOG.warn("There are {} entries in the samples database for Date: {}, Flowcell: {}, Lane: {}, Index: {}".format(str(len(samples)),*key.split("_")))
        while len(samples) > 1:
            keep = []
            for sample in samples:
                LOG.info("Project: {}, sample: {}, yield (M): {}, _id: {}".format(sample.get("sample_prj"),
                                                                                  sample.get("barcode_name"),
                                                                                  sample.get("read_count"),
                                                                                  sample.get("_id")))
                if not query_yes_no("Do you want to include this sample in the report?"):
                    LOG.info("Excluding sample run with _id: {}".format(sample.get("_id")))
                else:
                    keep.append(sample)
            samples = keep
        
        sample_d[key] = samples
    
    # Populate the sample run list with the verified sample runs
    sample_dict = [s[0] for s in sample_d.values()]
        
    rst_dict = dict()
    rst_dict.update(fc_dict[0])
    rst_dict.update(prj_dict)
    rst_dict["pdffile"] = "{}_{}_flowcell_summary.pdf".format(rst_dict["project_name"], rst_dict["FC_name"])
    rst_dict["rstfile"] = "{}.rst".format(os.path.splitext(rst_dict["pdffile"])[0])
    
    q30_cutoff = int(rst_dict["sample_q30_cutoff"])
    phix_cutoff = float(rst_dict["phix_cutoff"])
    lane_yield_cutoff = int(rst_dict["lane_yield_cutoff"])
    rst_dict["lane_yield_cutoff"] = _round_read_count_in_millions(lane_yield_cutoff)
    
    # Loop samples and build the sample information table
    snt_header = ['SciLifeLab ID','Submitted ID']
    syt_header = ['SciLifeLab ID','Lane','Barcode']
    syt_header.append('Read{}s (M)'.format(' pair' if rst_dict['is_paired'] else ''))
    sqt_header = ['SciLifeLab ID','Lane','Barcode','Avg Q','Q30 (%)','PhiX error rate (%)']
    lyt_header = ['Lane',
                  'Read{}s (M)'.format(' pair' if rst_dict['is_paired'] else '')]
    
    sample_name_table = [snt_header] + [[sample["scilifeid"], 
                                         sample["customerid"]] for sample in sample_dict]
    
    sample_yield_table = [syt_header] + [[sample["scilifeid"], 
                                          sample["lane"], 
                                          sample["barcode"],
                                          sample["read_count"]] for sample in sample_dict]
    
    sample_quality_table = [sqt_header] + [[sample["scilifeid"], 
                                            sample["lane"], 
                                            sample["barcode"],
                                            sample["avgQ"],
                                            "{}{}".format(sample["q30"],
                                                          "*" if float(sample["q30"]) < q30_cutoff else ""),
                                            "{}{}".format(sample["phix_error_rate"],
                                                          "*" if float(sample["phix_error_rate"]) > phix_cutoff else ""),
                                            ] for sample in sample_dict]
    
    lane_yield_table = [lyt_header] + [[lane,
                                        "{}{}".format(_round_read_count_in_millions(yld),
                                                      "*" if prj_dict["application"] == "Finished library" and yld < lane_yield_cutoff else "")] for lane,yld in fc_dict[0]["lane_yields"].items()]
    
    # Sort the tables by smaple and lane
    snt = [sample_name_table[0]] 
    for n in sorted(sample_name_table[1:], key=operator.itemgetter(0,1)): 
        if n not in snt:
            snt.append(n)
    sample_name_table = snt
    sample_yield_table = [sample_yield_table[0]] + sorted(sample_yield_table[1:], key=operator.itemgetter(0,2,3))
    lane_yield_table = [lane_yield_table[0]] + sorted(lane_yield_table[1:], key=operator.itemgetter(0))
    sample_quality_table = [sample_quality_table[0]] + sorted(sample_quality_table[1:], key=operator.itemgetter(0,2,3))

    # Write final output to reportlab and rst files
    output_data["debug"].write(json.dumps({'s_param': [rst_dict], 'sample_runs':{s["name"]:s["barcode_name"] for s in sample_dict}}))

    make_rest_note(rst_dict["rstfile"], 
                   tables={'name': sample_name_table, 'sample_yield': sample_yield_table, 'lane_yield': lane_yield_table, 'quality': sample_quality_table}, 
                   report="sample_report", **rst_dict)
    
    return output_data

def _get_qc_thresholds(params, config):
    """Get the specified QC thresholds from the config
    """
    info = {}
    
    # Get the PhiX error rate cutoff
    info['phix_cutoff'] = float(config.get("qc","phix_error_rate_threshold"))
    
    # Get the expected lane yield
    [hiseq_ho, hiseq_rm, miseq] = [False,False,False]
    lane_yield = None
    instr = params.get("instrument_version")
    if instr == "MiSeq":
        miseq = True
    elif instr.startswith("HiSeq"):
        if params.get("run_mode") == "RapidRun":
            hiseq_rm = True
        else:
            hiseq_ho = True
            
    if hiseq_ho:
        lane_yield = int(config.get("qc","hiseq_ho_lane_yield"))
    elif hiseq_rm:
        lane_yield = int(config.get("qc","hiseq_rm_lane_yield"))
    elif miseq:
        lane_yield = int(config.get("qc","miseq_lane_yield"))
    info['lane_yield_cutoff'] = lane_yield

    # Get the sample quality value cutoff
    cycles = params.get("num_cycles")
    pctq30 = 0
    for level in [250,150,100,50]:
        if cycles >= level:
            if miseq:
                pctq30 = int(config.get("qc","miseq_q30_{}".format(str(level))))
            elif hiseq_ho:
                pctq30 = int(config.get("qc","hiseq_ho_q30_{}".format(str(level))))
            elif hiseq_rm:
                pctq30 = int(config.get("qc","hiseq_rm_q30_{}".format(str(level))))
            break
    info['sample_q30_cutoff'] = pctq30
    
    return info

def _get_flowcell_info(fc_con, fc, project_name=None):
    info = {}
    info["FC_name"] = fc
    info["FC_id"] = fc_con.get_run_info(fc).get("Flowcell")
    info["FC_position"] = fc_con.get_run_parameters(fc).get("FCPosition")
    info["start_date"] = fc_con.get_start_date(fc)
    # Get instrument
    info['instrument_version'] = fc_con.get_instrument_type(fc)
    info['instrument_id'] = fc_con.get_instrument(fc)
    # Get run mode
    info["run_mode"] = fc_con.get_run_mode(fc)
    info["is_paired"] = fc_con.is_paired_end(fc)
    if info["is_paired"] is None:
        LOG.warn("Could not determine run setup for flowcell {}. Will assume paired-end.".format(fc))
        info["is_paired"] = True
    info["is_dual_index"] = fc_con.is_dual_index(fc)
    info["clustered"] = fc_con.get_clustered(fc)
    info["run_setup"] = fc_con.get_run_setup(fc)
    info["num_cycles"] = fc_con.num_cycles(fc)
    info["lane_yields"] = fc_con.get_lane_yields(fc,project_name)
    info["phix_error_rate"] = {lane: fc_con.get_phix_error_rate(fc,lane) for lane in info["lane_yields"].keys()}
    info["software"] = fc_con.get_demultiplex_software(fc)
    
    return info
    
    
def _get_project_info(project, pcon):
    info = {}
    info["project_name"] = project.get("project_name")
    info["customer_reference"] = project.get("customer_reference")
    info["application"] = project.get("application")
    info["lanes_ordered"] = project.get("details",{}).get("sequence_units_ordered_(lanes)")
    info["no_samples"] = project.get("no_of_samples")
    info["uppnex_project_id"] = project.get("uppnex_id")
    info["project_id"] = project.get("project_id")
    info["open_date"] = project.get("open_date")
    info["samples"] = project.get("samples")
    info["m_ordered"] = pcon.get_ordered_amount(info["project_name"],samples=info["samples"])
    
    return info

def _get_project_sample_info(project):
    """Get the sample information stored in the project database
    """
    samples = []
    for sample in project.get("samples",{}).values():
        LOG.debug("got sample {}".format(sample.get("scilife_name")))
        sample_dict = {}
        sample_dict['scilifeid'] = sample.get('scilife_name')
        sample_dict['customerid'] = sample.get('customer_name')
        sample_dict['m_ordered'] = sample.get('reads_requested_(millions)',project.get('min_m_reads_per_sample_ordered'))
        s = sample.get('incoming_QC_status','')
        if s is not None:
            s = (s.lower() == 'p' or s.lower() == 'passed' or s.lower() == 'true' or s.lower() == 't')
        sample_dict['incoming_qc_status'] = s
        
        for item in ['start_date','finish_date']:
            try:
                s = datetime.datetime.strptime(sample.get('initial_qc_',{}).get(item),'%Y-%m-%d')
            except:
                s = None
            sample_dict['incoming_qc_{}'.format(item)] = s
            
        sample_dict['overall_status'] = sample.get('status')
        sample_dict['m_reads_sequenced'] = sample.get('m_reads_sequenced')
        sample_dict['sample_type'] = sample.get('details',{}).get('sample_type')
        sample_dict['fragment_size_customer'] = sample.get('details',{}).get('customer_average_fragment_length')
        
        # Get the library prep information
        preps = []
        for lbl, prep in sample.get('library_prep',{}).items():
            prep_dict = {}
            prep_dict['label'] = lbl
            s = prep.get('prep_status')
            if s is not None:
                s = (s.lower() == 'p' or s.lower() == 'passed' or s.lower() == 'true' or s.lower() == 't')
            prep_dict['prep_qc_status'] = s
            
            for item in ['start_date','finished_date']:
                try:
                    s = datetime.datetime.strptime(prep.get('prep_{}'.format(item)),'%Y-%m-%d')
                except:
                    s = None
                prep_dict['prep_{}'.format(item.replace('finished','finish'))] = s
        
            prep_dict['fragment_size'] = prep.get('average_size_bp')
            prep_dict['reagent_label'] = ";".join(prep.get('reagent_labels',[]))
            
            # Get the library validation information
            libvals = prep.get('library_validation',{})
            prep_dict['fragment_size'] = libvals.get('average_size_bp',prep_dict['fragment_size'])
            for val in libvals.values():
                if not type(val) is dict:
                    continue
                for item in ['start_date','finish_date']:
                    key = 'prep_qc_{}'.format(item)
                    try:
                        s = datetime.datetime.strptime(val[item],'%Y-%m-%d')
                        if  key not in prep_dict or prep_dict[key] < s:
                            prep_dict[key] = s
                    except:
                        prep_dict[key] = prep_dict.get(key)
                        
                prep_dict['fragment_size'] = val.get('average_size_bp',prep_dict['fragment_size'])
            
            # Get the sequencing run information
            runs = []
            for lbl, run_metric in prep.get('sample_run_metrics',{}).items():
                run_dict = {}
                run_dict['sample_run_name'] = lbl
                if type(run_metric) is str:
                    run_dict['sample_run_id'] = run_metric
                else:
                    run_dict['sample_run_id'] = run_metric.get('sample_run_metrics_id')
                    for key, name in [('sequencing_finish_date','run_finish_date'),
                                      ('dillution_and_pooling_start_date','pool_start_date'),
                                      ('sequencing_run_QC_finished','run_qc_finish_date'),]:
                        try:
                            s = datetime.datetime.strptime(run_metric.get(key),'%Y-%m-%d')
                        except:
                            s = None
                        run_dict[name] = s
                runs.append(run_dict)
            prep_dict['sequencing_runs'] = runs
            preps.append(prep_dict)
        sample_dict['library_preps'] = preps
        samples.append(sample_dict)
    return samples

def _get_sample_info(sample_run_list, project, fc_con, **kw):
    
    samples = []
    # Loop samples and build the sample information 
    for sample in sample_run_list: 
        LOG.debug("working on sample '{}', sample run metrics name '{}', id '{}'".format(sample.get("barcode_name"), 
                                                                                         sample.get("name"), 
                                                                                         sample.get("_id")))
        # Get the project sample name corresponding to the sample run
        project_sample = project.get("samples",{}).get(sample.get("project_sample_name"))
        if project_sample:
            # FIXME: Is this really necessary? There doesn't seem to be any consequence if the ids don't match
            LOG.debug("project sample run metrics mapping found: '{}' : '{}'".format(sample["name"], project_sample["scilife_name"]))
            # Set project_sample_d: a dictionary mapping from sample run metrics name to sample run metrics database id
            project_sample_d = _set_project_sample_dict(project_sample, project.get("source","gdocs"))
            if not project_sample_d:
                LOG.warn("No sample_run_metrics information for sample '{}', barcode name '{}', id '{}'\n\tProject summary information {}".format(sample["name"], 
                                                                                                                                                  sample["barcode_name"], 
                                                                                                                                                  sample["_id"], 
                                                                                                                                                  project_sample))
            # Check if sample run metrics name present in project database: if so, verify that database ids are consistent
            if sample["name"] not in project_sample_d.keys():
                LOG.warn("no such sample run metrics '{}' in project sample run metrics dictionary".format(sample["name"]) )
            else:
                if sample["_id"] == project_sample_d[sample["name"]]:
                    LOG.debug("project sample run metrics mapping found: '{}' : '{}'".format(sample["name"], 
                                                                                             project_sample_d[sample["name"]]))
                else:
                    LOG.warn("inconsistent mapping for '{}': '{}' != '{}' (project summary id)".format(sample["name"], 
                                                                                                       sample["_id"], 
                                                                                                       project_sample_d[sample["name"]]))
            sample['customer_name'] = project_sample.get("customer_name")
            
        # No project sample found. Manual upload to database necessary.
        else:
            sample['customer_name'] = None
            LOG.warn("No project sample name found for sample run name '{}'".format(sample["barcode_name"]))
            LOG.info("Please run 'pm qc upload-qc FLOWCELL_ID --extensive-matching' to update project sample names ")
            LOG.info("or 'pm qc update --sample_prj PROJECT_NAME --names BARCODE_TO_SAMPLE_MAP to update project sample names.")
            LOG.info("Please refer to the pm documentation for examples.")
            query_ok(force=kw.get("force", False))
        
        # Get read counts, possible overridden on the command line
        if kw.get("bc_count"):
            read_count = _round_read_count_in_millions(_get_bc_count(sample["barcode_name"], _literal_eval_option(kw.get("bc_count")), sample))
        else:
            read_count = _round_read_count_in_millions(sample.get("bc_count"))
        
        # Get quality score from demultiplex stats, if that fails
        # (which it shouldn't), fall back on fastqc data.
        fc = "{}_{}".format(sample["date"],
                            sample["flowcell"])
        (avg_quality_score, pct_q30_bases) = fc_con.get_barcode_lane_statistics(project.get("project_name"), 
                                                                                sample.get("barcode_name"), 
                                                                                fc, 
                                                                                sample["lane"])
        if not avg_quality_score:
            avg_quality_score = calc_avg_qv(sample) 
        if not avg_quality_score:
            LOG.warn("Setting average quality failed for sample {}, id {}".format(sample.get("name"), 
                                                                                  sample.get("_id")))
        if not pct_q30_bases:
            LOG.warn("Setting % of >= Q30 Bases (PF) failed for sample {}, id {}".format(sample.get("name"), 
                                                                                         sample.get("_id")))
        
        # Get phix error rate, possible overridden on the command line
        if kw.get("phix"):
            phix_rate = _get_phix_error_rate(sample["lane"], _literal_eval_option(kw.get("phix")))
        else:
            phix_rate = fc_con.get_phix_error_rate(fc,sample["lane"]) 
        
        sample_item = dict()
        sample_item["_id"] = sample.get("_id")
        sample_item["name"] = sample.get("name")
        sample_item["scilifeid"] = sample.get("project_sample_name")
        sample_item["sample_prj"] = sample.get("sample_prj")
        sample_item["customerid"] = sample.get("customer_name")
        sample_item["barcode_name"] = sample.get("barcode_name")
        sample_item["lane"] = sample.get("lane")
        sample_item["barcode"] = sample.get("sequence")
        sample_item["date"] = sample.get("date")
        sample_item["flowcell"] = sample.get("flowcell")
        sample_item["read_count"] = read_count
        sample_item["avgQ"] = avg_quality_score
        sample_item["q30"] = pct_q30_bases
        sample_item["phix_error_rate"] = phix_rate
        samples.append(sample_item)
    
    return samples

def _exclude_sample_id(exclude_sample_ids, sample_name, barcode_seq):
    """Check whether we should exclude a sample id.

    :param exclude_sample_ids: dictionary of sample:barcode pairs
    :param sample_name: project sample name
    :param barcode_seq: the barcode sequence

    :returns: True if exclude, False otherwise
    """
    if exclude_sample_ids and sample_name in exclude_sample_ids.keys():
        if exclude_sample_ids[sample_name]:
            if barcode_seq in exclude_sample_ids[sample_name]:
                LOG.info("excluding sample '{}' with barcode '{}' from project report".format(sample_name, barcode_seq))
                return True
            else:
                LOG.info("keeping sample '{}' with barcode '{}' in sequence report".format(sample_name, barcode_seq))
                return False
        else:
            LOG.info("excluding sample '{}' from project report".format(sample_name))
            return True


def _set_sample_table_values(sample_name, project_sample, barcode_seq, ordered_million_reads, param):
    """Set the values for a sample that is to appear in the final table.

    :param sample_name: string identifier of sample
    :param project_sample: project sample dictionary from project summary database
    :param barcode_seq: barcode sequence
    :param ordered_million_reads: the number of ordered reads
    :param param: project parameters

    :returns: vals, a dictionary of table values
    """
    prjs_to_table = {'ScilifeID':'scilife_name', 'SubmittedID':'customer_name', 'MSequenced':'m_reads_sequenced'}#, 'MOrdered':'min_m_reads_per_sample_ordered', 'Status':'status'}
    vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
    # Set status
    vals['Status'] = project_sample.get("status", "N/A")
    if ordered_million_reads:
        param["ordered_amount"] = _get_ordered_million_reads(sample_name, ordered_million_reads)
    vals['MOrdered'] = param["ordered_amount"]
    vals['BarcodeSeq'] = barcode_seq
    vals.update({k:"N/A" for k in vals.keys() if vals[k] is None or vals[k] == ""})
    return vals

def data_delivery_note(**kw):
    """Create an easily parseable information file with information about the data delivery
    """
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    
    project_name = kw.get('project_name',None)
    flowcell = kw.get('flowcell',None)
    LOG.debug("Generating data delivery note for project {}{}.".format(project_name,' and flowcell {}'.format(flowcell if flowcell else '')))
    
    # Get a connection to the project and sample databases
    p_con = ProjectSummaryConnection(**kw)
    assert p_con, "Could not connect to project database"
    s_con = SampleRunMetricsConnection(**kw)
    assert s_con, "Could not connect to sample database"
    
    # Get the entry for the project and samples from the database
    LOG.debug("Fetching samples from sample database")
    samples = s_con.get_samples(sample_prj=project_name, fc_id=flowcell)
    LOG.debug("Got {} samples from database".format(len(samples)))
    
    # Get the customer sample names from the project database
    LOG.debug("Fetching samples from project database")
    project_samples = p_con.get_entry(project_name, "samples")
    customer_names = {sample_name:sample.get('customer_name','N/A') for sample_name, sample in project_samples.items()}
    
    data = [['SciLifeLab ID','Submitted ID','Flowcell','Lane','Barcode','Read','Path','MD5','Size (bytes)','Timestamp']]
    for sample in samples:
        sname = sample.get('project_sample_name','N/A')
        cname = customer_names.get(sname,'N/A')
        fc = sample.get('flowcell','N/A')
        lane = sample.get('lane','N/A')
        barcode = sample.get('sequence','N/A')
        if 'raw_data_delivery' not in sample:
            data.append([sname,cname,'','','','','','','',''])
            continue
        delivery = sample['raw_data_delivery']
        tstamp = delivery.get('timestamp','N/A')
        for read, file in delivery.get('files',{}).items():
            data.append([sname,
                         cname,
                         fc,
                         lane,
                         barcode,
                         read,
                         file.get('path','N/A'),
                         file.get('md5','N/A'),
                         file.get('size_in_bytes','N/A'),
                         tstamp,])
    
    # Write the data to a csv file
    outfile = "{}{}_data_delivery.csv".format(project_name,'_{}'.format(flowcell) if flowcell else '')
    LOG.debug("Writing delivery data to {}".format(outfile))
    with open(outfile,"w") as outh:
        csvw = csv.writer(outh)
        for row in data:
            csvw.writerow(row)
    
    # Write Texttable formatted output to stdout
    tt = texttable.Texttable(180)
    tt.add_rows(data)
    output_data['stdout'].write(tt.draw())
        
    return output_data
    

def project_status_note(project_name=None, username=None, password=None, url=None,
                        use_ps_map=True, use_bc_map=False, check_consistency=False,
                        ordered_million_reads=None, uppnex_id=None, customer_reference=None,
                        exclude_sample_ids={}, project_alias=None, sample_aliases={},
                        projectdb="projects", samplesdb="samples", flowcelldb="flowcells",
                        include_all_samples=False, flat_table=False, **kw):
    """Make a project status note. Used keywords:

    :param project_name: project name
    :param user: db user name
    :param password: db password
    :param url: db url
    :param use_ps_map: use project summary mapping
    :param use_bc_map: use project to barcode name mapping
    :param check_consistency: check consistency between mappings
    :param ordered_million_reads: number of ordered reads in millions
    :param uppnex_id: the uppnex id
    :param customer_reference: customer project name
    :param exclude_sample_ids: exclude some sample ids from project note
    :param project_alias: project alias name
    :param sample_aliases: sample alias names
    :param projectdb: project db name
    :param samplesdb: samples db name
    :param flowcelldb: flowcells db name
    :param include_all_samples: include all samples in report
    :param flat_table: Just create a simple tab-separated version of the table instead of the fancy pdf
    """

    kw["project_name"] = project_name
    kw["username"] = username
    kw["password"] = password
    kw["url"] = url
    kw["uppnex_id"] = uppnex_id
    kw["customer_reference"] = customer_reference
    kw["project_alias"] = project_alias
    kw["projectdb"] = projectdb
    kw["samplesdb"] = samplesdb
    kw["flowcelldb"] = flowcelldb
    kw["use_ps_map"] = use_ps_map
    kw["use_bc_map"] = use_bc_map
    kw["check_consistency"] = check_consistency
    kw["ordered_million_reads"] = ordered_million_reads
    kw["exclude_sample_ids"] = exclude_sample_ids
    kw["sample_aliases"] = sample_aliases
    kw["include_all_samples"] = include_all_samples
    kw["flat_table"] = flat_table
    
    # parameters
    parameters = {
        "project_name" : project_name,
        "finished" : "Not finished, or cannot yet assess if finished.",
        }

    output_data, tables, param = _project_status_note_table(**kw)

    if not flat_table:
        # Set report paragraphs
        #paragraphs = project_note_paragraphs()
        #headers = project_note_headers()

        #paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
        #make_note("{}_project_summary.pdf".format(project_name), headers, paragraphs, **param)
        make_rest_note("{}_project_summary.rst".format(param["project_name"]), 
                       report="project_report", 
                       tables=tables, 
                       **param)

    else:
        # Write tab-separated output
        sample_table[0].insert(0,'ProjectID')
        table_cols = [sample_table[0].index(col) for col in ['ProjectID', 'ScilifeID', 'SubmittedID', 'BarcodeSeq', 'MSequenced']]
        outfile = "{}_project_summary.csv".format(project_name)
        with open(outfile,"w") as outh:
            csvw = csv.writer(outh)
            for i,sample in enumerate(sample_table):
                if i > 0:
                    sample.insert(0,project_name)
                data = [str(sample[col]) for col in table_cols]
                csvw.writerow(data)
                output_data['stdout'].write("{}\n".format("\t".join(data)))

    param.update({k:"N/A" for k in param.keys() if param[k] is None or param[k] ==  ""})
    output_data["debug"].write(json.dumps({'param':param, 'sample_table': tables["sample_table"], 'flowcell_table': tables["flowcell_table"]}))

    return output_data


def _project_status_note_table(**kw):

    # mapping project_summary to parameter keys
    ps_to_parameter = {"scilife_name":"scilife_name", "customer_name":"customer_name", "project_name":"project_name"}
    # mapping project sample to table
    table_keys = ['ScilifeID', 'SubmittedID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status']

    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    """
    # Connect and run
    s_con = SampleRunMetricsConnection(dbname=samplesdb, username=username, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(dbname=flowcelldb, username=username, password=password, url=url)
    p_con = ProjectSummaryConnection(dbname=projectdb, username=username, password=password, url=url)

    #Get the information source for this project
    source = p_con.get_info_source(project_name)

    # Get project summary from project database
    sample_aliases = _literal_eval_option(sample_aliases, default={})
    prj_summary = p_con.get_entry(project_name)
    if not prj_summary:
        LOG.warn("No such project '{}'".format(project_name))
        return
    LOG.debug("Working on project '{}'.".format(project_name))
    
    # Get sample run list and loop samples to make mapping sample -> {sampleruns}
    sample_run_list = _set_sample_run_list(project_name, flowcell=None, project_alias=project_alias, s_con=s_con)
    samples = {}
    for s in sample_run_list:
        prj_sample = p_con.get_project_sample(project_name, s.get("project_sample_name", None))
        if prj_sample:
            sample_name = prj_sample['project_sample'].get("scilife_name", None)
            s_d = {s["name"] : {'sample':sample_name, 'id':s["_id"]}}
            samples.update(s_d)
        else:
            if s["barcode_name"] in sample_aliases:
                s_d = {sample_aliases[s["barcode_name"]] : {'sample':sample_aliases[s["barcode_name"]], 'id':s["_id"]}}
                samples.update(s_d)
            else:
                s_d = {s["name"]:{'sample':s["name"], 'id':s["_id"], 'barcode_name':s["barcode_name"]}}
                LOG.warn("No mapping found for sample run:\n  '{}'".format(s_d))

    # Convert to mapping from desired sample name to list of aliases
    # Less important for the moment; one solution is to update the
    # Google docs summary table to use the P names
    sample_dict = prj_summary['samples']
    param.update({key:prj_summary.get(ps_to_parameter[key], None) for key in ps_to_parameter.keys()})
    param["ordered_amount"] = param.get("ordered_amount", p_con.get_ordered_amount(project_name, samples=sample_dict))
    param['customer_reference'] = param.get('customer_reference', prj_summary.get('customer_reference'))
    param['uppnex_project_id'] = param.get('uppnex_project_id', prj_summary.get('uppnex_id'))

    # Override database values if options passed at command line
    if uppnex_id:
        param["uppnex_project_id"] = uppnex_id
    if customer_reference:
        param["customer_reference"] = customer_reference
    """
    
    data = _collect_status_note_data(**kw)
    if len(data) == 0:
        return output_data
    
    prj_dict, fc_dict, sample_dict = data[0:3]
   
    # Process options
    prj_dict["m_ordered"] = _literal_eval_option(kw.get("ordered_million_reads"), default=prj_dict.get("m_ordered"))
    exclude_sample_ids = _literal_eval_option(kw.get("exclude_sample_ids"), default={})

    # Create a list of lanes containing this project's samples on each flowcell
    fc_lanes = {}
    for sample in sample_dict:
        key = "{}_{}".format(sample["date"],sample["flowcell"])
        if key not in fc_lanes:
            fc_lanes[key] = []
        fc_lanes[key] = list(set(fc_lanes[key] + [sample["lane"]])) 

    # Create the flowcell table
    fct_header = ["Start date", "Flowcell position", "Flowcell id", "Lane", "Yield (M)", "PhiX error rate (%)"]
    fc_table = [fct_header] + [[fc["start_date"],
                                fc["FC_position"],
                                fc["FC_id"],
                                lane,
                                "{}{}".format(_round_read_count_in_millions(fc["lane_yields"][lane]),
                                              "*" if fc["lane_yields"][lane] < fc["lane_yield_cutoff"] else ""),
                                "{}{}".format(fc["phix_error_rate"][lane],
                                              "*" if fc["phix_error_rate"][lane] > fc["phix_cutoff"] else "")] for fc in fc_dict for lane in fc_lanes[fc["FC_name"]]]
    fc_table = [fc_table[0]] + sorted(fc_table[1:], key=operator.itemgetter(0,1,2,3))
    
    # Create the sample table
    st_header = ["SciLifeLab ID", "Submitted ID", "Total yield (M)", "Status"]
    sample_table = [st_header] + [[sample.get("scilife_name"),
                                   sample.get("customer_name"),
                                   sample.get("m_reads_sequenced","0"),
                                   sample.get("status","N/A")] for sample in prj_dict.get("samples",{}).values()]
    sample_table = [sample_table[0]] + sorted(sample_table[1:], key=operator.itemgetter(0,1,3,2))
    
    rst_dict = prj_dict
    rst_dict["flowcell_dicts"] = fc_dict
    rst_dict["sample_dicts"] = sample_dict
    
    """ 
    ## Start collecting the data
    sample_table = []
    samples_excluded = []
    all_passed = True
    last_library_preps = p_con.get_latest_library_prep(project_name)
    last_library_preps_srm = [x for l in last_library_preps.values() for x in l]
    LOG.debug("Looping through sample map that maps project sample names to sample run metrics ids")
    for k,v in samples.items():
        LOG.debug("project sample '{}' maps to '{}'".format(k, v))
        if not include_all_samples:
            if v['sample'] not in last_library_preps.keys():
                LOG.info("No library prep information for sample {}; keeping in report".format(v['sample']))
            else:
                if k not in last_library_preps_srm:
                    LOG.info("Sample run {} ('{}') is not latest library prep ({}) for project sample {}: excluding from report".format(k, v["id"], ",".join(list(set(last_library_preps[v['sample']].values()))), v['sample']))
                    continue
        else:
            pass

        if re.search("Unexpected", k):
            continue
        barcode_seq = s_con.get_entry(k, "sequence")
        # Exclude sample id?
        if _exclude_sample_id(exclude_sample_ids, v['sample'], barcode_seq):
            samples_excluded.append(v['sample'])
            continue
        # Get the project sample name from the sample run and set table values
        project_sample = sample_dict[v['sample']]
        vals = _set_sample_table_values(v['sample'], project_sample, barcode_seq, ordered_million_reads, param)
        if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
        sample_table.append([vals[k] for k in table_keys])

    # Loop through samples in sample_dict for which there is no sample run information
    samples_in_table_or_excluded = list(set([x[0] for x in sample_table])) + samples_excluded
    samples_not_in_table = list(set(sample_dict.keys()) - set(samples_in_table_or_excluded))
    for sample in samples_not_in_table:
        if re.search("Unexpected", sample):
            continue
        project_sample = sample_dict[sample]
        # Set project_sample_d: a dictionary mapping from sample run metrics name to sample run metrics database id
        project_sample_d = _set_project_sample_dict(project_sample, source)
        if project_sample_d:
            for k,v in project_sample_d.iteritems():
                barcode_seq = s_con.get_entry(k, "sequence")
                vals = _set_sample_table_values(sample, project_sample, barcode_seq, ordered_million_reads, param)
                if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
                sample_table.append([vals[k] for k in table_keys])
        else:
            barcode_seq = None
            vals = _set_sample_table_values(sample, project_sample, barcode_seq, ordered_million_reads, param)
            if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
            sample_table.append([vals[k] for k in table_keys])
    if all_passed: param["finished"] = 'Project finished.'
    sample_table.sort()
    sample_table = list(sample_table for sample_table,_ in itertools.groupby(sample_table))
    sample_table.insert(0, ['ScilifeID', 'SubmittedID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
    """
    return output_data, {"flowcell_table": fc_table, "sample_table": sample_table}, rst_dict

