"""Module delivery_notes - code for generating delivery reports and notes"""
import re
import itertools
from cStringIO import StringIO

from scilifelab.db.statusdb import SampleRunMetricsConnection, ProjectSummaryConnection, FlowcellRunMetricsConnection
from scilifelab.report import sequencing_success
from scilifelab.report.rl import make_note, concatenate_notes, sample_note_paragraphs, sample_note_headers, project_note_paragraphs, project_note_headers, make_sample_table
import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)


def sample_status_note(project_id=None, flowcell_id=None, user=None, password=None, url=None,
                       use_ps_map=True, use_bc_map=False, check_consistency=False, 
                       ordered_million_reads=None, uppnex_id=None, customer_reference=None,
                       qcinfo=False, **kw):
    """Make a sample status note. Used keywords:

    :param project_id: project id
    :param flowcell_id: flowcell id
    :param user: db user name
    :param password: db password
    :param url: db url
    :param use_ps_map: use project summary mapping
    :param use_bc_map: use project to barcode name mapping
    :param check_consistency: check consistency between mappings
    :param ordered_million_reads: number of ordered reads in millions
    :param uppnex_id: the uppnex id
    :param customer_reference: customer project name
    :param qcinfo: flag to print qc info
    """
    ## Cutoffs
    cutoffs = {
        "phix_err_cutoff" : 2.0,
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
    ## Write qcinfo if needed
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    if qcinfo:
        qcdata["stdout"].write("*** Quality stats ***\n")
        qcdata["stdout"].write("Scilifelab ID\tPhiXError\tAvgQV\n")
        
    ## Connect and run
    s_con = SampleRunMetricsConnection(username=user, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(username=user, password=password, url=url)
    p_con = ProjectSummaryConnection(username=user, password=password, url=url)
    paragraphs = sample_note_paragraphs()
    headers = sample_note_headers()
    project = p_con.get_entry(project_id)
    notes = []
    if not project:
        LOG.warn("No such project '{}'".format(project_id))
        return
    samples = p_con.map_srm_to_name(project_id, include_all=False, fc_id=flowcell_id, use_ps_map=use_ps_map, use_bc_map=use_bc_map, check_consistency=check_consistency)
    for k,v  in samples.items():
        s_param = {}
        LOG.debug("working on sample '{}', sample run metrics name '{}', id '{}'".format(v["sample"], k, v["id"]))
        s_param.update(parameters)
        if not v['id'] is None:
            if not s_con.name_fc_view[k].value == flowcell_id:
                LOG.debug("skipping sample '{}' since it isn't run on flowcell {}".format(k, flowcell_id))
                continue
        else:
            if re.search("NOSRM", k):
                LOG.warn("No sample run metrics information for project sample '{}'".format(k.strip("NOSRM_")))
                continue
        s = s_con.get_entry(k)
        s_param.update({key:s[srm_to_parameter[key]] for key in srm_to_parameter.keys()})
        fc = "{}_{}".format(s["date"], s["flowcell"])
        s_param["phix_error_rate"] = fc_con.get_phix_error_rate(str(fc), s["lane"])
        s_param['avg_quality_score'] = s_con.calc_avg_qv(s["name"])
        if qcinfo:
            self.app._output_data["stdout"].write("{}\t{}\t{}\n".format(s["barcode_name"], s_param["phix_error_rate"], s_param["avg_quality_score"]))
        s_param['rounded_read_count'] = round(float(s_param['rounded_read_count'])/1e6,1) if s_param['rounded_read_count'] else None
        s_param['ordered_amount'] = s_param.get('ordered_amount', p_con.get_ordered_amount(project_id))
        s_param['customer_reference'] = s_param.get('customer_reference', project.get('customer_reference'))
        s_param['uppnex_project_id'] = s_param.get('uppnex_project_id', project.get('uppnex_id'))
        if ordered_million_reads:
            s_param["ordered_amount"] = ordered_million_reads
        if uppnex_id:
            s_param["uppnex_project_id"] = uppnex_id
        if customer_reference:
            s_param["customer_reference"] = customer_reference
        s_param['customer_name'] = project['samples'].get(v["sample"], {}).get("customer_name", None)
        s_param['success'] = sequencing_success(s_param, cutoffs)
        s_param.update({k:"N/A" for k in s_param.keys() if s_param[k] is None or s_param[k] ==  ""})
        notes.append(make_note("{}_{}_{}.pdf".format(s["barcode_name"], s["date"], s["flowcell"]), headers, paragraphs, **s_param))
    concatenate_notes(notes, "{}_{}_{}_sample_summary.pdf".format(project_id, s["date"], s["flowcell"]))
    return output_data

def project_status_note(project_id=None, user=None, password=None, url=None,
                        use_ps_map=True, use_bc_map=False, check_consistency=False, 
                        ordered_million_reads=None, uppnex_id=None, customer_reference=None,
                        exclude_sample_ids=[], **kw):
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
        
    ## Connect and run
    s_con = SampleRunMetricsConnection(username=user, password=password, url=url)
    fc_con = FlowcellRunMetricsConnection(username=user, password=password, url=url)
    p_con = ProjectSummaryConnection(username=user, password=password, url=url)
    paragraphs = project_note_paragraphs()
    headers = project_note_headers()
    param = parameters
    project = p_con.get_entry(project_id)
    if not project:
        LOG.warn("No such project '{}'".format(project_id))
        return
    LOG.debug("Working on project '{}'.".format(project_id))
    samples = p_con.map_srm_to_name(project_id, use_ps_map=use_ps_map, use_bc_map=use_bc_map, check_consistency=check_consistency)
    sample_list = project['samples']
    param.update({key:project.get(ps_to_parameter[key], None) for key in ps_to_parameter.keys()})
    param["ordered_amount"] = param.get("ordered_amount", p_con.get_ordered_amount(project_id))
    param['customer_reference'] = param.get('customer_reference', project.get('customer_reference'))
    param['uppnex_project_id'] = param.get('uppnex_project_id', project.get('uppnex_id'))
    if ordered_million_reads:
        param["ordered_amount"] = ordered_million_reads
    if uppnex_id:
        param["uppnex_project_id"] = uppnex_id
    if customer_reference:
        param["customer_reference"] = customer_reference
    if not param["ordered_amount"]:
        param["ordered_amount"] = ordered_million_reads
    ## Start collecting the data
    sample_table = []
    all_passed = True
    LOG.debug("Looping through sample map that maps project sample names to sample run metrics ids")
    for k,v in samples.items():
        LOG.debug("project sample '{}' maps to '{}'".format(k, v))
        if re.search("Unexpected", k):
            continue
        if exclude_sample_ids and v['sample'] in exclude_sample_ids[0].split():
            LOG.info("excluding sample '{}' from project report".format(v['sample']))
            continue
        project_sample = sample_list[v['sample']]
        vals = {x:project_sample.get(prjs_to_table[x], None) for x in prjs_to_table.keys()}
        ## Set status
        vals['Status'] = project_sample.get("status", "N/A")
        vals['MOrdered'] = param["ordered_amount"]
        vals['BarcodeSeq'] = s_con.get_entry(k, "sequence")
        vals.update({k:"N/A" for k in vals.keys() if vals[k] is None or vals[k] == ""})
        if vals['Status']=="N/A" or vals['Status']=="NP": all_passed = False
        sample_table.append([vals[k] for k in table_keys])
    if all_passed: param["finished"] = 'Project finished.'
    sample_table.sort()
    sample_table = list(sample_table for sample_table,_ in itertools.groupby(sample_table))
    sample_table.insert(0, ['ScilifeID', 'CustomerID', 'BarcodeSeq', 'MSequenced', 'MOrdered', 'Status'])
    paragraphs["Samples"]["tpl"] = make_sample_table(sample_table)
    make_note("{}_project_summary.pdf".format(project_id), headers, paragraphs, **param)

