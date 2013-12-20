import copy
import datetime

from scilifelab.db.statusdb import ProjectSummaryConnection, SampleRunMetricsConnection, FlowcellRunMetricsConnection

def report(ctrl, project_name, **kw):
    
    import ipdb
    ipdb.set_trace()
    data = fetch_project_data(ctrl, project_name, kw.get('sample_name'), kw.get('flowcell_id'), **kw)
    import pprint
    pprint.pprint(data)

def fetch_project_data(ctrl, project_name, sample_name=None, flowcell_id=None, **kw):
    """Fetch all relevant information from StatusDB for a project in order to make a sequencing report.
    If sample_name or flowcell_id are specified, the data will be pruned to those sample(s) or flowcell(s).
    Returns a data structure representing the information
    """
    
    pcon = ProjectSummaryConnection(**kw)
    scon = SampleRunMetricsConnection(**kw)
    fcon = FlowcellRunMetricsConnection(**kw)
    
    for obj, db in [(pcon,"project"),(scon,"sample"),(fcon,"flowcell")]:
        assert obj, "Could not connect to {} database in StatusDB".format(db)
    
    # Wrap the arguments in lists in order to iterate
    if type(project_name) is not list:
        project_name = [project_name]
    if sample_name and type(sample_name) is not list:
        sample_name = [sample_name]
    if flowcell_id and type(flowcell_id) is not list:
        flowcell_id = [flowcell_id]
    
    data = {}
    for pname in project_name:
        
        project = pcon.get_entry(pname)
        if not project:
            ctrl.warn("No such project '{}'".format(pname))
            continue
        if project.get('source') != 'lims':
            ctrl.warn("The source for data for project {} is not LIMS. This type of report can only be generated for data from LIMS. Please refer to older versions of pm".format(pname))
            continue
        
        # Fetch the information frpom the project database
        pdata = _get_project_info(project,sample_name,flowcell_id)
        
        # Fetch sample run information from the samples database
        sdata = _get_sample_run_info(scon,pdata.get("project_id"),sample_name,flowcell_id)
        fcs = ["_".join([sr.get("date").strftime("%y%m%d"),sr.get("flowcell")]) for sruns in sdata.values() for sr in sruns]
        
        # Fetch run information from the flowcell database for the flowcells listed among the sample runs
        for fcid in list(set(fcs)):
            if flowcell_id and fcid in flowcell_id:
                continue
            flowcell, sample_results = _get_flowcell_run_info(fcon.get_entry(fcid),pname,sample_name)
            if "flowcells" not in pdata:
                pdata["flowcells"] = {}
            pdata["flowcells"][fcid] = flowcell
            
            # Add the sample results under the corresponding sample_run object
            for result in sample_results.values():
                key = "_".join([result.get("lane"),fcid,result.get("Index")])
                for sname, sruns in sdata.items():
                    for i, srun in enumerate(sruns):
                        if srun.get("name") != key:
                            continue
                        # Verify that the project sample name and the sample id match before joining the results
                        if not result.get("sampleid").startswith(srun.get("project_sample_name")):
                            ctrl.error("Flowcell results for {} does not have a sample name that match the project_sample_name " \
                                       "of the corresponding sample run ({} vs {})".format(key,
                                                                                           result.get("sampleid"),
                                                                                           srun.get("project_sample_name")))
                            continue
                        sdata[sname][i].update(result)
                        key = None
                        break
                    if not key:
                        break
                # If key is not None, no match was found for the sample among the sample runs
                if key:
                    ctrl.warn("Flowcell results were found for sample {} ({}) " \
                              "but no corresponding sample run was found".format(result.get("sampleid"),
                                                                                 key))
        
        # Add the sample run object to the corresponding sample in the project
        for sname, sruns in sdata.items():
            if sname not in pdata.get("samples",{}):
                ctrl.warn("Sample run for sample '{}' found in samples database but sample is missing from project database".format(sname))
                pdata["samples"][sname] = {}
            pdata["samples"][sname]["sample_runs"] = sruns
            
        data[pname] = pdata
        
    return data

def _get_flowcell_lane_info(demux_stats, project, sample_name=None):
    """Summarize the demultiplex info related to a project
    """
    lanes = {}
    samples = {}
    # Loop over the entries in Demultiplex_Stats and store information under lane keys
    for entry in demux_stats:
        lane = entry.get("Lane")
        index = entry.get("Index")
        if lane not in lanes:
            lanes[lane] = {"reads": 0, "projects": [], "samples": []}
            
        # Sum the reads across the lane
        reads = _parse_int(entry.get("# Reads"))
        lanes[lane]["reads"] += reads
        if index == "Undetermined":
            lanes[lane]["undetermined"] = reads
        
        # Replace any '__' in project name with '.'
        pname = entry.get("Project","").replace("__",".")
        lanes[lane]["projects"] = list(set(lanes[lane]["projects"] + [pname]))
        
        # Skip any entries not belonging to the project we're interested in
        if pname != project:
            continue
        
        # Skip any samples that are in the list of samples to ignore
        sname = entry.get("Sample ID","").replace("__",".")
        if sample_name and sname in sample_name:
            continue
        lanes[lane]["samples"] = list(set(lanes[lane]["samples"] + [sname]))
        
        # Store the sample details in the sample dict
        
        skey = "_".join([lane,index])
        samples[skey] = {"sampleid": sname,
                         "lane": lane, 
                         "reads": reads, 
                         "Index": index,
                         "q30": float(entry.get("% of >= Q30 Bases (PF)","-1").replace(",",".")),
                         "avgq": float(entry.get("Mean Quality Score (PF)","-1").replace(",",".")),
                         "pct0mm": float(entry.get("% Perfect Index Reads","-1").replace(",",".")),
                         "pct1mm": float(entry.get("% One Mismatch Reads (Index)","-1").replace(",",".")),
                         "yield": _parse_int(entry.get("Yield (Mbases)"))*1e6
                         }
    return lanes, samples

def _get_flowcell_run_details(runinfo=None, runparameters=None):
    """Get the flowcell run details from the RunInfo and RunParameters
    dictionaries
    """
    info = {}
    if runinfo:
        info.update({f:runinfo.get(f) for f in ["Number",
                                                "Instrument",
                                                "Id",
                                                "Flowcell"]})
        info["Reads"] = {r.get("Number"): {"NumCycles":r.get("NumCycles"), 
                                           "IsIndexedRead":_translate_QC_flag("IsIndexedRead")} for r in runinfo.get("Reads",[])}
        info["Date"] = _parse_date(runinfo.get("Date"),"%y%m%d")
        
    if runparameters:
        info.update({f:runparameters.get(f) for f in ["RunMode",
                                                      "RTAVersion",
                                                      "ClusteringChoice",
                                                      "ApplicationName",
                                                      "ApplicationVersion",
                                                      "FCPosition"]})
        info["RunStartDate"] = _parse_date(runparameters.get("RunStartDate"),"%y%m%d")
        info["PairEndFC"] = _translate_QC_flag(runparameters.get("PairEndFC"))
                                                      
    return info

def _get_flowcell_run_info(flowcell, project, sample_name=None):
    """Get the flowcell run information from the flowcell database
    """
    info = {f:flowcell.get(f) for f in ["name",
                                        "run_setup"]}
    lanes, samples = _get_flowcell_lane_info(flowcell.get("illumina",{}).get("Demultiplex_Stats",{}).get("Barcode_lane_statistics",[]),project)
    info["lanes"] = lanes
    info["runinfo"] = _get_flowcell_run_details(flowcell.get("RunInfo"),flowcell.get("RunParameters"))
    info["runsummary"] = _get_flowcell_run_summary(flowcell.get("illumina",{}).get("run_summary"))
    
    return info, samples 

def _get_flowcell_run_summary(summary):
    """Get the run summary of the flowcell that was fetched from LIMS
    """
    info = {}
    if not summary:
        return info

    for lane, data in summary.items():
        lane_info = {f:data.get(f) for f in ["Cluster Density (K/mm^2) R1",
                                             "Cluster Density (K/mm^2) R2",
                                             "% Bases >=Q30 R1",
                                             "% Bases >=Q30 R2",
                                             "%PF R1",
                                             "%PF R2",
                                             "Avg Q Score R1",
                                             "Avg Q Score R2",
                                             "Clusters PF R1",
                                             "Clusters PF R2",
                                             "Yield PF (Gb) R1",
                                             "Yield PF (Gb) R2",
                                             "Clusters Raw R1",
                                             "Clusters Raw R2",
                                             "% Error Rate R1",
                                             "% Error Rate R2"]}
        lane_info["qc"] = _translate_QC_flag(summary.get("qc"))
        info.update({lane: lane_info})
        
    return info
                    
def _get_library_prep_info(prep, fcid=None):
    """Get the library prep information stored in the project database
    """ 
    info = {f:_translate_QC_flag(prep.get(f)) for f in ["prep_status"]}
    info.update({f:_parse_date(prep.get(f)) for f in ["prep_start_date",
                                                      "prep_finished_date"]})
    info["reagent_labels"] = ";".join(prep.get("reagent_labels",[]))
    info["library_validation"] = [_get_library_validation_info(v) for v in prep.get("library_validation",{}).values()]
    info["sample_run_metrics"] = {k:_get_project_sample_run_metrics_info(v) for k,v in prep.get("sample_run_metrics",{}).items() if _use_flowcell(k,fcid)}
    return info
    
def _get_library_validation_info(val):
    """Get the library validation information stored in the project database
    """
    info = {f:_parse_date(val.get(f)) for f in ["start_date",
                                                "finish_date"]}
    info.update({f:val.get(f) for f in ["average_size_bp"]})
    return info

def _get_project_info(project, sample_name=None, flowcell_id=None):
    info = {f:project.get(f) for f in ["project_name",
                                       "customer_reference",
                                       "application",
                                       "no_of_samples",
                                       "uppnex_id",
                                       "project_id"]}
    info.update({f:_parse_date(project.get(f)) for f in ["open_date",
                                                         "sequencing_finished"]})
    
    details = project.get('details')
    info.update({f:details.get(f) for f in ["best_practice_bioinformatics",
                                            "bioinformatic_qc",
                                            "custom_capture_design_id",
                                            "customer_project_reference",
                                            "portal_id",
                                            "sequence_units_ordered_(lanes)",
                                            "type",]})
    info.update({f:_parse_date(details.get(f)) for f in ["queued"]})
    info["samples"] = {k:_get_project_sample_info(v,flowcell_id) for k,v in project.get("samples",{}).items() if sample_name is None or k not in sample_name}
    return info
     
def _get_project_sample_info(sample, fcid=None):
    """Get the sample information stored in the project database
    """
        
    # Get general sample information
    info = {f:sample.get(f) for f in ["scilife_name",
                                      "customer_name",
                                      "m_reads",
                                      "reads_requested_(millions)",
                                      "status"]}
    # Check the alternative key for sequenced amount
    if info["m_reads"] is None:
        info["m_reads"] = sample.get("m_reads_sequenced")
        
    info.update({f:_translate_QC_flag(sample.get(f)) for f in ["incoming_QC_status",
                                                               "status"]})
    info.update({f:_parse_date(sample.get(f)) for f in ["incoming_qc_start_date",
                                                        "incoming_qc_finish_date"]})
    details = sample.get("details",{})
    info.update({f:details.get(f) for f in ["sample_type"]})
    
    # Get the library prep information
    info["library_prep"] = {k:_get_library_prep_info(v,fcid) for k,v in sample.get("library_prep",{}).items() if type(v) is dict}
    
    return info
   
def _get_project_sample_run_metrics_info(run):
    """Get the sample run metrics from the project database
    """
    info = {f:_parse_date(run.get(f)) for f in ["dillution_and_pooling_start_date",
                                                "sequencing_run_QC_finished",
                                                "sequencing_start_date",
                                                "sequencing_finish_date"]}
    info.update({f:run.get(f) for f in ["sample_run_metrics_id"]})
    return info

def _get_sample_run_info(scon, project_id, sample_name=None, flowcell_id=None):
    """Get the sample run information for samples belonging to a project
    """
    samples = scon.get_samples_by_project_id(project_id)
    sinfo = {}
    for sample in samples:
        if sample_name and sample.get("project_sample_name") in sample_name:
            continue
        if not _use_flowcell(sample.get("name"),flowcell_id):
            continue
        info = {f:sample.get(f) for f in ["name",
                                          "lane",
                                          "flowcell",
                                          "project_id",
                                          "sequence",
                                          "project_sample_name",
                                          "_id"]}
        info["date"] = _parse_date(sample.get("date"),"%y%m%d")
        if not info["project_sample_name"] in sinfo:
            sinfo[info["project_sample_name"]] = []
        sinfo[info["project_sample_name"]].append(info)
        
    return sinfo

def _parse_date(val,fmt="%Y-%m-%d"):
    try:
        val = datetime.datetime.strptime(val,fmt)
    except TypeError, ValueError:
        val = None
    return val

def _parse_int(val):
    """Parse a int value with comma separated thousands
    """
    return 0 if not val else int(val.replace(",",""))

def _translate_QC_flag(val):
    if val is not None and type(val) == str:
        val = (val.lower() == 'p' or val.lower() == 'passed' or val.lower() == 'true' or val.lower() == 't' or val.lower() == 'y')
    return val

def _use_flowcell(s,haystack):
    if haystack:
        # Join date and flowcell id
        needle = "_".join(s.split("_")[1:3])
        for straw in haystack:
            s = straw.split("_")
            if (len(s) > 1 and needle == "_".join([s[0],s[-1]])) or \
                (len(s) == 1 and needle.endswith(s[0])):
                return False
    return True
