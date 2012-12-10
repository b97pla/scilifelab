"""Module delivery_notes - code for generating delivery reports and notes"""
import os
import re
import itertools
import ast
import json
import math
from cStringIO import StringIO
from collections import Counter
from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection, calc_avg_qv
from scilifelab.report import sequencing_success
from scilifelab.report.sphinx import make_rest_note
from scilifelab.report.rl import make_note, concatenate_notes, sample_note_paragraphs, sample_note_headers, project_note_paragraphs, project_note_headers, make_sample_table
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

# http://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
def _round_read_count_in_millions(n):
    """Round absolute read counts to million reads"""
    LOG.debug("Rounding read count: got {}".format(n))
    if n is None:
        return None
    if n == 0:
        return 0
    round_factor = [2,2,1]
    millidx = max(0, min(len(round_factor) - 1, int(math.floor(math.log10(abs(n))/3.0))))
    return round(float(n)/10**(6),round_factor[millidx])

def _get_ordered_million_reads(sample_name, ordered_million_reads):
    """Retrieve ordered million reads for sample

    :param sample_name: sample name
    :param ordered_million_reads: parsed option passed to application

    :returns: ordered number of reads or None"""
    if isinstance(ordered_million_reads, dict):
        if sample_name in ordered_million_reads:
            return ordered_million_reads[sample_name]
        else:
            return ordered_million_reads.get("default", -1)
    else:
        return ordered_million_reads

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
    """Assert name of flowcell: "[A-Z0-9]+XX"

    :param flowcell: flowcell id
    
    :returns: boolean
    """
    if flowcell is None:
        return True
    if not re.match("[A-Z0-9]+XX$", flowcell):
        return False
    return True

def sample_status_note(project_id=None, flowcell=None, username=None, password=None, url=None,
                       ordered_million_reads=None, uppnex_id=None, customer_reference=None, bc_count=None,
                       project_alias=[], projectdb="projects", samplesdb="samples", flowcelldb="flowcells",
                       phix=None, **kw):
    """Make a sample status note. Used keywords:

    :param project_id: project id
    :param flowcell: flowcell id
    :param username: db username
    :param password: db password
    :param url: db url
    :param ordered_million_reads: number of ordered reads in millions
    :param uppnex_id: the uppnex id
    :param customer_reference: customer project name
    :param project_alias: project alias name
    :param phix: phix error rate
    """
    ## Cutoffs
    cutoffs = {
        "phix_err_cutoff" : 2.0,
        "qv_cutoff" : 30,
        }
    
    ## parameters
    parameters = {
        "project_name" : None,
        "start_date" : None,
        "FC_id" : None,
        "scilifelab_name" : None,
        "rounded_read_count" : None,
        "phix_error_rate" : None,
        "avg_quality_score" : None,
        "success" : None,
        }
        ## key mapping from sample_run_metrics to parameter keys
    srm_to_parameter = {"project_name":"sample_prj", "FC_id":"flowcell", 
                        "scilifelab_name":"barcode_name", "start_date":"date", "rounded_read_count":"bc_count"}
    
    LOG.debug("got parameters {}".format(parameters))
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    if not _assert_flowcell_format(flowcell):
        LOG.warn("Wrong flowcell format {}; skipping. Please use the flowcell id (format \"[A-Z0-9]+XX\")".format(flowcell) )
        return output_data

    output_data["stdout"].write("\nQuality stats\n")
    output_data["stdout"].write("************************\n")
    output_data["stdout"].write("PhiX error cutoff: > {:3}\n".format(cutoffs['phix_err_cutoff']))
    output_data["stdout"].write("QV cutoff        : < {:3}\n".format(cutoffs['qv_cutoff']))
    output_data["stdout"].write("************************\n\n")
    output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format("Scilifelab ID", "Lane", "PhiXError", "ErrorStatus", "AvgQV", "QVStatus"))
    output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format("=============", "====", "=========", "===========", "=====", "========"))
    # Connect and run
    s_con = SampleRunMetricsConnection(dbname=samplesdb, username=username, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(dbname=flowcelldb, username=username, password=password, url=url)
    p_con = ProjectSummaryConnection(dbname=projectdb, username=username, password=password, url=url)
    paragraphs = sample_note_paragraphs()
    headers = sample_note_headers()
    project = p_con.get_entry(project_id)
    notes = []
    s_param_out = {}
    if not project:
        LOG.warn("No such project '{}'".format(project_id))
        return output_data
    sample_run_list = s_con.get_samples(sample_prj=project_id, fc_id=flowcell)
    if project_alias:
        project_alias = ast.literal_eval(project_alias)
        for p_alias in project_alias:
            sample_run_list_tmp = s_con.get_samples(sample_prj=p_alias, fc_id=flowcell)
            if sample_run_list_tmp:
                sample_run_list.extend(sample_run_list_tmp)
    if ordered_million_reads:
        if os.path.exists(ordered_million_reads):
            with open(ordered_million_reads) as fh:
                ordered_million_reads = json.load(fh)
        else:
            ordered_million_reads = ast.literal_eval(ordered_million_reads)
    if bc_count:
        if os.path.exists(bc_count):
            with open(bc_count) as fh:
                bc_count = json.load(fh)
        else:
            bc_count = ast.literal_eval(bc_count)
    if len(sample_run_list) == 0:
        LOG.warn("No samples for project '{}', flowcell '{}'. Maybe there are no sample run metrics in statusdb?".format(project_id, flowcell))
        return output_data
    sample_count = Counter([x.get("barcode_name") for x in sample_run_list])
    for s in sample_run_list:
        s_param = {}
        LOG.debug("working on sample '{}', sample run metrics name '{}', id '{}'".format(s.get("barcode_name", None), s.get("name", None), s.get("_id", None)))
        s_param.update(parameters)
        s_param.update({key:s[srm_to_parameter[key]] for key in srm_to_parameter.keys()})
        fc = "{}_{}".format(s.get("date"), s.get("flowcell"))
        s_param["phix_error_rate"] = phix if phix else fc_con.get_phix_error_rate(str(fc), s["lane"])
        s_param['avg_quality_score'] = calc_avg_qv(s)
        if not s_param['avg_quality_score']:
            LOG.warn("Calculation of average quality failed for sample {}, id {}".format(s.get("name"), s.get("_id")))
        err_stat = "OK"
        qv_stat = "OK"
        if s_param["phix_error_rate"] > cutoffs["phix_err_cutoff"]:
            err_stat = "HIGH"
        elif s_param["phix_error_rate"] == -1:
            err_stat = "N/A"
        if s_param["avg_quality_score"] < cutoffs["qv_cutoff"]:
            qv_stat = "LOW"
        output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format(s["barcode_name"], s["lane"], s_param["phix_error_rate"], err_stat, s_param["avg_quality_score"], qv_stat))
        s_param['ordered_amount'] = s_param.get('ordered_amount', p_con.get_ordered_amount(project_id))
        s_param['customer_reference'] = s_param.get('customer_reference', project.get('customer_reference'))
        s_param['uppnex_project_id'] = s_param.get('uppnex_project_id', project.get('uppnex_id'))
        if ordered_million_reads:
            s_param["ordered_amount"] = _get_ordered_million_reads(s["barcode_name"], ordered_million_reads)
        if bc_count:
            s_param["rounded_read_count"] = _round_read_count_in_millions(_get_bc_count(s["barcode_name"], bc_count, s))
        else:
            s_param["rounded_read_count"] = _round_read_count_in_millions(s_param["rounded_read_count"])
        if uppnex_id:
            s_param["uppnex_project_id"] = uppnex_id
        if customer_reference:
            s_param["customer_reference"] = customer_reference
        ## FIX ME: This is where we need a key in SampleRunMetrics that provides a mapping to a project sample name
        project_sample = p_con.get_project_sample(project_id, s["barcode_name"])
        if project_sample:
            project_sample_item = project_sample['project_sample']
            if "library_prep" in project_sample_item.keys():
                project_sample_d = {x:y for d in [v["sample_run_metrics"] for k,v in project_sample_item["library_prep"].iteritems()] for x,y in d.iteritems()}
            else:
                project_sample_d = {x:y for x,y in project_sample_item.get("sample_run_metrics", {}).iteritems()}
                if not project_sample_item.get("sample_run_metrics", {}):
                    LOG.warn("No sample_run_metrics information for sample '{}', barcode name '{}', id '{}'\n\tProject summary information {}".format(s["name"], s["barcode_name"], s["_id"], project_sample_item))
            if s["name"] not in project_sample_d.keys():
                LOG.warn("'{}' not found in project sample run metrics for project".format(s["name"]) )
            else:
                if s["_id"] == project_sample_d[s["name"]]:
                    LOG.debug("project sample run metrics mapping found: '{}' : '{}'".format(s["name"], project_sample_d[s["name"]]))
                else:
                    LOG.warn("inconsistent mapping for '{}': '{}' != '{}' (project summary id)".format(s["name"], s["_id"], project_sample_d[s["name"]]))
            s_param['customer_name'] = project_sample_item.get("customer_name", None)
        else:
            s_param['customer_name'] = None
            LOG.warn("No project sample name found for sample run name '{}'".format(s["barcode_name"]))
        s_param['success'] = sequencing_success(s_param, cutoffs)
        s_param.update({k:"N/A" for k in s_param.keys() if s_param[k] is None or s_param[k] ==  "" or s_param[k] == -1.0})
        if sample_count[s.get("barcode_name")] > 1:
            outfile = "{}_{}_{}_{}.pdf".format(s["barcode_name"], s["date"], s["flowcell"], s["lane"])
        else:
            outfile = "{}_{}_{}.pdf".format(s["barcode_name"], s["date"], s["flowcell"])
        notes.append(make_note(outfile, headers, paragraphs, **s_param))
        make_rest_note(outfile.replace(".pdf", ".rst"), **s_param)
        s_param_out[s_param["scilifelab_name"]] = s_param
    output_data["debug"].write(json.dumps({'s_param': s_param_out, 'sample_runs':{s["name"]:s["barcode_name"] for s in sample_run_list}}))
    concatenate_notes(notes, "{}_{}_{}_sample_summary.pdf".format(project_id, s.get("date", None), s.get("flowcell", None)))
    return output_data

def project_status_note(project_id=None, username=None, password=None, url=None,
                        use_ps_map=True, use_bc_map=False, check_consistency=False,
                        ordered_million_reads=None, uppnex_id=None, customer_reference=None,
                        exclude_sample_ids={}, project_alias=None, sample_aliases={},
                       projectdb="projects", samplesdb="samples", flowcelldb="flowcells", **kw):
    """Make a project status note. Used keywords:

    :param project_id: project id
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
    """
    ## parameters
    parameters = {
        "project_name" : project_id,
        "finished" : "Not finished, or cannot yet assess if finished.",
        }
    ## mapping project_summary to parameter keys
    ps_to_parameter = {"scilife_name":"scilife_name", "customer_name":"customer_name", "project_name":"project_id"}
    ## mapping project sample to table
    table_keys = ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status']
    prjs_to_table = {'ScilifeID':'scilife_name', 'CustomerID':'customer_name', 'MSequenced':'m_reads_sequenced'}#, 'MOrdered':'min_m_reads_per_sample_ordered', 'Status':'status'}

    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    ## Connect and run
    s_con = SampleRunMetricsConnection(dbname=samplesdb, username=username, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(dbname=flowcelldb, username=username, password=password, url=url)
    p_con = ProjectSummaryConnection(dbname=projectdb, username=username, password=password, url=url)
    paragraphs = project_note_paragraphs()
    headers = project_note_headers()
    param = parameters
    prj_summary = p_con.get_entry(project_id)
    if sample_aliases:
        if os.path.exists(sample_aliases):
            with open(sample_aliases) as fh:
                sample_aliases = json.load(fh)
        else:
            sample_aliases = ast.literal_eval(sample_aliases)
    if not prj_summary:
        LOG.warn("No such project '{}'".format(project_id))
        return
    LOG.debug("Working on project '{}'.".format(project_id))
    sample_run_list = s_con.get_samples(sample_prj=project_id)
    # If project aliases, extend sample_run_list with results
    if project_alias:
        project_alias = ast.literal_eval(project_alias)
        for p_alias in project_alias:
            sample_run_list_tmp = s_con.get_samples(sample_prj=p_alias)
            if sample_run_list_tmp:
                sample_run_list.extend(sample_run_list_tmp)
    samples = {}
    for s in sample_run_list:
        prj_sample = p_con.get_project_sample(project_id, s["barcode_name"])
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
    ## Convert to mapping from desired sample name to list of aliases
    ## Less important for the moment; one solution is to update the
    ## Google docs summary table to use the P names
    sample_list = prj_summary['samples']
    param.update({key:prj_summary.get(ps_to_parameter[key], None) for key in ps_to_parameter.keys()})
    param["ordered_amount"] = param.get("ordered_amount", p_con.get_ordered_amount(project_id))
    param['customer_reference'] = param.get('customer_reference', prj_summary.get('customer_reference'))
    param['uppnex_project_id'] = param.get('uppnex_project_id', prj_summary.get('uppnex_id'))
    if ordered_million_reads:
        if os.path.exists(ordered_million_reads):
            with open(ordered_million_reads) as fh:
                ordered_million_reads = json.load(fh)
        else:
            ordered_million_reads = ast.literal_eval(ordered_million_reads)
    if uppnex_id:
        param["uppnex_project_id"] = uppnex_id
    if customer_reference:
        param["customer_reference"] = customer_reference
    if exclude_sample_ids:
        if os.path.exists(exclude_sample_ids):
            with open(exclude_sample_ids) as fh:
                exclude_sample_ids = json.load(fh)
        else:
            exclude_sample_ids = ast.literal_eval(exclude_sample_ids)
    ## Start collecting the data
    sample_table = []
    all_passed = True
    LOG.debug("Looping through sample map that maps project sample names to sample run metrics ids")
    for k,v in samples.items():
        LOG.debug("project sample '{}' maps to '{}'".format(k, v))
        if re.search("Unexpected", k):
            continue
        barcode_seq = s_con.get_entry(k, "sequence")
        if exclude_sample_ids and v['sample'] in exclude_sample_ids.keys():
            if exclude_sample_ids[v['sample']]:
                if barcode_seq in exclude_sample_ids[v['sample']]:
                    LOG.info("excluding sample '{}' with barcode '{}' from project report".format(v['sample'], barcode_seq))
                    continue
                else:
                    LOG.info("keeping sample '{}' with barcode '{}' in sequence report".format(v['sample'], barcode_seq))
            else:
                LOG.info("excluding sample '{}' from project report".format(v['sample']))
                continue
        project_sample = sample_list[v['sample']]
        vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
        ## Set status
        vals['Status'] = project_sample.get("status", "N/A")
        if ordered_million_reads:
            param["ordered_amount"] = _get_ordered_million_reads(v['sample'], ordered_million_reads)
        vals['MOrdered'] = param["ordered_amount"]
        vals['BarcodeSeq'] = barcode_seq
        vals.update({k:"N/A" for k in vals.keys() if vals[k] is None or vals[k] == ""})
        if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
        sample_table.append([vals[k] for k in table_keys])
    if all_passed: param["finished"] = 'Project finished.'
    sample_table.sort()
    sample_table = list(sample_table for sample_table,_ in itertools.groupby(sample_table))
    sample_table.insert(0, ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
    paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
    make_note("{}_project_summary.pdf".format(project_id), headers, paragraphs, **param)
    make_rest_note("{}_project_summary.rst".format(project_id), **param)
    param.update({k:"N/A" for k in param.keys() if param[k] is None or param[k] ==  ""})
    output_data["debug"].write(json.dumps({'param':param, 'table':sample_table}))
    return output_data

