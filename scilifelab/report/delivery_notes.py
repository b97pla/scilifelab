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
from scilifelab.utils.misc import query_ok
from scilifelab.report import sequencing_success
from scilifelab.report.rst import make_sample_rest_notes, make_rest_note
from scilifelab.report.rl import make_note, concatenate_notes, sample_note_paragraphs, sample_note_headers, project_note_paragraphs, project_note_headers, make_sample_table
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

# Instrument configurations for insertion into report
instrument = {
    'SN1025': {
        'instrument_alias':'Smeagol',
        'instrument' : 'HiSeq 2000',
        },
    'sN188': {
        'instrument_alias':'Eowyn',# ?
        'instrument' : 'HiSeq 2000',
        },
    'SN7001298': {
        'instrument_alias':'Deagol',
        'instrument' : 'HiSeq 2000',
        },
    'SN7001301': {
        'instrument_alias':'Eomer',
        'instrument' : 'HiSeq 2500',
        },
    'SN7001362': {
        'instrument_alias':'Aragorn',
        'instrument' : 'HiSeq 2500',
        },
    'default': {
        'instrument_alias':'',
        'instrument' : 'NN',
        },
    'M00275' : {
        'instrument_alias':'Gloin',
        'instrument' : 'MiSeq 2500',
        }
    }
# Software versions used in data production. Instrument specific?
software_versions = {
    'baseconversion_version' : 'OLB v1.9',
    'casava_version' : 'CASAVA v1.8.2'
    }


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
    """Assert name of flowcell: "[A-Z0-9]+XX"

    :param flowcell: flowcell id
    
    :returns: boolean
    """
    if flowcell is None:
        # Can this really be right?!?
        return True
    if not re.match("[A-Z0-9]+XX$", flowcell):
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

def _set_project_sample_dict(project_sample_item):
    """Set a project sample dict, mapping a project sample to sample run metrics if present in project summary.

    :param project_sample_item: a project sample item

    :returns: project_sample_d or empty dict
    """
    project_sample_d = {}
    if "library_prep" in project_sample_item.keys():
        project_sample_d = {x:y for d in [v.get("sample_run_metrics", {}) for k,v in project_sample_item["library_prep"].iteritems()] for x,y in d.iteritems()}
    else:
        project_sample_d = {x:y for x,y in project_sample_item.get("sample_run_metrics", {}).iteritems()}
        if not project_sample_item.get("sample_run_metrics", {}):
            LOG.warn("No sample_run_metrics information for sample '{}'".format(project_sample_item))
    return project_sample_d


def sample_status_note(project_name=None, flowcell=None, username=None, password=None, url=None,
                       ordered_million_reads=None, uppnex_id=None, customer_reference=None, bc_count=None,
                       project_alias=[], projectdb="projects", samplesdb="samples", flowcelldb="flowcells",
                       phix=None, **kw):
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
    """
    # Cutoffs
    cutoffs = {
        "phix_err_cutoff" : 2.0,
        "qv_cutoff" : 30,
        }
    
    # parameters
    parameters = {
        "project_name" : None,
        "start_date" : None,
        "FC_id" : None,
        "scilifelab_name" : None,
        "rounded_read_count" : None,
        "phix_error_rate" : None,
        "avg_quality_score" : None,
        "pct_q30_bases" : None,
        "success" : None,
        "run_mode":None,
        }
    # key mapping from sample_run_metrics to parameter keys
    srm_to_parameter = {"project_name":"sample_prj", "FC_id":"flowcell", 
                        "scilifelab_name":"barcode_name", "start_date":"date", "rounded_read_count":"bc_count"}
    
    LOG.debug("got parameters {}".format(parameters))
    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    if not _assert_flowcell_format(flowcell):
        LOG.warn("Wrong flowcell format {}; skipping. Please use the flowcell id (format \"[A-Z0-9]+XX\")".format(flowcell) )
        return output_data
    output_data = _update_sample_output_data(output_data, cutoffs)

    # Connect and run
    s_con = SampleRunMetricsConnection(dbname=samplesdb, username=username, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(dbname=flowcelldb, username=username, password=password, url=url)
    p_con = ProjectSummaryConnection(dbname=projectdb, username=username, password=password, url=url)

    # Set up paragraphs
    paragraphs = sample_note_paragraphs()
    headers = sample_note_headers()

    # Get project
    project = p_con.get_entry(project_name)
    if not project:
        LOG.warn("No such project '{}'".format(project_name))
        return output_data

    # Set samples list
    sample_run_list = _set_sample_run_list(project_name, flowcell, project_alias, s_con)
    if len(sample_run_list) == 0:
        LOG.warn("No samples for project '{}', flowcell '{}'. Maybe there are no sample run metrics in statusdb?".format(project_name, flowcell))
        return output_data
    
    # Set options
    ordered_million_reads = _literal_eval_option(ordered_million_reads)
    bc_count = _literal_eval_option(bc_count)
    phix = _literal_eval_option(phix)

    # Count number of times a sample has been run on a flowcell; if several, make lane-specific reports
    sample_count = Counter([x.get("barcode_name") for x in sample_run_list])

    # Loop samples and collect information
    s_param_out = []
    for s in sample_run_list:
        s_param = {}
        LOG.debug("working on sample '{}', sample run metrics name '{}', id '{}'".format(s.get("barcode_name", None), s.get("name", None), s.get("_id", None)))
        s_param.update(parameters)
        s_param.update({key:s[srm_to_parameter[key]] for key in srm_to_parameter.keys()})
        fc = "{}_{}".format(s.get("date"), s.get("flowcell"))
        # Get instrument
        try:
            s_param.update(instrument[fc_con.get_instrument(str(fc))])
        except:
            LOG.warn("Failed to set instrument and software versions for flowcell {} in report due to missing RunInfo -> Instrument field in statusdb. Either rerun 'pm qc update-qc' or search-and-replace 'NN' in the sample report.".format(fc))
            s_param.update(instrument['default'])
        # Get run mode
        s_param["run_mode"] = fc_con.get_run_mode(str(fc))
        s_param.update(software_versions)
        s_param["phix_error_rate"] = fc_con.get_phix_error_rate(str(fc), s["lane"])
        if phix:
            s_param["phix_error_rate"] = _get_phix_error_rate(s["lane"], phix)
        # Get quality score from demultiplex stats, if that fails
        # (which it shouldn't), fall back on fastqc data.
        (avg_quality_score, pct_q30_bases) = fc_con.get_barcode_lane_statistics(project_name, s.get("barcode_name"), fc, s["lane"])
        s_param['avg_quality_score'] = avg_quality_score if avg_quality_score else calc_avg_qv(s) 
        if not s_param['avg_quality_score']:
            LOG.warn("Setting average quality failed for sample {}, id {}".format(s.get("name"), s.get("_id")))
        s_param['pct_q30_bases'] = pct_q30_bases
        if not s_param['pct_q30_bases']:
            LOG.warn("Setting % of >= Q30 Bases (PF) failed for sample {}, id {}".format(s.get("name"), s.get("_id")))
        # Compare phix error and qv to cutoffs
        err_stat = "OK"
        qv_stat = "OK"
        if s_param["phix_error_rate"] > cutoffs["phix_err_cutoff"]:
            err_stat = "HIGH"
        elif s_param["phix_error_rate"] == -1:
            err_stat = "N/A"
        if s_param["avg_quality_score"] < cutoffs["qv_cutoff"]:
            qv_stat = "LOW"
        output_data["stdout"].write("{:>18}\t{:>6}\t{:>12}\t{:>12}\t{:>12}\t{:>12}\n".format(s["barcode_name"], s["lane"], s_param["phix_error_rate"], err_stat, s_param["avg_quality_score"], qv_stat))

        # Update/set remaning sample run parameters, falling back on project defaults if *key* is missing
        s_param['ordered_amount'] = s_param.get('ordered_amount', p_con.get_ordered_amount(project_name))
        s_param['customer_reference'] = s_param.get('customer_reference', project.get('customer_reference'))
        s_param['uppnex_project_id'] = s_param.get('uppnex_project_id', project.get('uppnex_id'))

        # Override database settings if options passed at command line
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

        # Get the project sample name corresponding to the sample run
        project_sample = p_con.get_project_sample(project_name, s.get("project_sample_name", None))
        if project_sample:
            LOG.debug("project sample run metrics mapping found: '{}' : '{}'".format(s["name"], project_sample["sample_name"]))
            project_sample_item = project_sample['project_sample']
            # Set project_sample_d: a dictionary mapping from sample run metrics name to sample run metrics database id
            project_sample_d = _set_project_sample_dict(project_sample_item)
            if not project_sample_d:
                LOG.warn("No sample_run_metrics information for sample '{}', barcode name '{}', id '{}'\n\tProject summary information {}".format(s["name"], s["barcode_name"], s["_id"], project_sample))
            # Check if sample run metrics name present in project database: if so, verify that database ids are consistent
            if s["name"] not in project_sample_d.keys():
                LOG.warn("no such sample run metrics '{}' in project sample run metrics dictionary".format(s["name"]) )
            else:
                if s["_id"] == project_sample_d[s["name"]]:
                    LOG.debug("project sample run metrics mapping found: '{}' : '{}'".format(s["name"], project_sample_d[s["name"]]))
                else:
                    LOG.warn("inconsistent mapping for '{}': '{}' != '{}' (project summary id)".format(s["name"], s["_id"], project_sample_d[s["name"]]))
            s_param['customer_name'] = project_sample_item.get("customer_name", None)

        # No project sample found. Manual upload to database necessary.
        else:
            s_param['customer_name'] = None
            LOG.warn("No project sample name found for sample run name '{}'".format(s["barcode_name"]))
            LOG.info("Please run 'pm qc upload-qc FLOWCELL_ID --extensive-matching' to update project sample names ")
            LOG.info("or 'pm qc update --sample_prj PROJECT_NAME --names BARCODE_TO_SAMPLE_MAP to update project sample names.")
            LOG.info("Please refer to the pm documentation for examples.")
            query_ok(force=kw.get("force", False))

        # Finally assess sequencing success, update parameters and set outputs
        s_param['success'] = sequencing_success(s_param, cutoffs)
        s_param.update({k:"N/A" for k in s_param.keys() if s_param[k] is None or s_param[k] ==  "" or s_param[k] == -1.0})
        if sample_count[s.get("barcode_name")] > 1:
            outfile = "{}_{}_{}_{}.pdf".format(s["barcode_name"], s["date"], s["flowcell"], s["lane"])
        else:
            outfile = "{}_{}_{}.pdf".format(s["barcode_name"], s["date"], s["flowcell"])
        s_param["outfile"] = outfile
        s_param_out.append(s_param)

    # Write final output to reportlab and rst files
    output_data["debug"].write(json.dumps({'s_param': s_param_out, 'sample_runs':{s["name"]:s["barcode_name"] for s in sample_run_list}}))
    notes = [make_note(headers=headers, paragraphs=paragraphs, **sp) for sp in s_param_out]
    rest_notes = make_sample_rest_notes("{}_{}_{}_sample_summary.rst".format(project_name, s.get("date", None), s.get("flowcell", None)), s_param_out)
    concatenate_notes(notes, "{}_{}_{}_sample_summary.pdf".format(project_name, s.get("date", None), s.get("flowcell", None)))
    return output_data


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
    prjs_to_table = {'ScilifeID':'scilife_name', 'CustomerID':'customer_name', 'MSequenced':'m_reads_sequenced'}#, 'MOrdered':'min_m_reads_per_sample_ordered', 'Status':'status'}
    vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
    # Set status
    vals['Status'] = project_sample.get("status", "N/A")
    if ordered_million_reads:
        param["ordered_amount"] = _get_ordered_million_reads(sample_name, ordered_million_reads)
    vals['MOrdered'] = param["ordered_amount"]
    vals['BarcodeSeq'] = barcode_seq
    vals.update({k:"N/A" for k in vals.keys() if vals[k] is None or vals[k] == ""})
    return vals

def project_status_note(project_name=None, username=None, password=None, url=None,
                        use_ps_map=True, use_bc_map=False, check_consistency=False,
                        ordered_million_reads=None, uppnex_id=None, customer_reference=None,
                        exclude_sample_ids={}, project_alias=None, sample_aliases={},
                        projectdb="projects", samplesdb="samples", flowcelldb="flowcells",
                        include_all_samples=False, **kw):
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
    """
    # parameters
    parameters = {
        "project_name" : project_name,
        "finished" : "Not finished, or cannot yet assess if finished.",
        }
    # mapping project_summary to parameter keys
    ps_to_parameter = {"scilife_name":"scilife_name", "customer_name":"customer_name", "project_name":"project_name"}
    # mapping project sample to table
    table_keys = ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status']

    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    # Connect and run
    s_con = SampleRunMetricsConnection(dbname=samplesdb, username=username, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(dbname=flowcelldb, username=username, password=password, url=url)
    p_con = ProjectSummaryConnection(dbname=projectdb, username=username, password=password, url=url)

    # Set report paragraphs
    paragraphs = project_note_paragraphs()
    headers = project_note_headers()
    # Set local param variable
    param = parameters
    
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
    param["ordered_amount"] = param.get("ordered_amount", p_con.get_ordered_amount(project_name))
    param['customer_reference'] = param.get('customer_reference', prj_summary.get('customer_reference'))
    param['uppnex_project_id'] = param.get('uppnex_project_id', prj_summary.get('uppnex_id'))

    # Override database values if options passed at command line
    if uppnex_id:
        param["uppnex_project_id"] = uppnex_id
    if customer_reference:
        param["customer_reference"] = customer_reference

    # Process options
    ordered_million_reads = _literal_eval_option(ordered_million_reads)
    exclude_sample_ids = _literal_eval_option(exclude_sample_ids, default={})

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
                    LOG.info("Sample run {} ('{}') is not latest library prep ({}) for project sample {}: excluding from report".format(k, v["id"], last_library_preps[v['sample']].values()[0], v['sample']))
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
        project_sample_d = _set_project_sample_dict(project_sample)
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
    sample_table.insert(0, ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
    paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
    make_note("{}_project_summary.pdf".format(project_name), headers, paragraphs, **param)
    make_rest_note("{}_project_summary.rst".format(project_name), sample_table=sample_table, report="project_report", **param)
    param.update({k:"N/A" for k in param.keys() if param[k] is None or param[k] ==  ""})
    output_data["debug"].write(json.dumps({'param':param, 'table':sample_table}))
    return output_data

