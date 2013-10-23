"""A module for writing lane QC data (typically from RTA) to google docs
"""

import os 
from bcbio.log import create_log_handler
from bcbio.pipeline.qcsummary import RTAQCMetrics
from bcbio.pipeline.flowcell import Flowcell
from bcbio.google import (_from_unicode,_to_unicode,get_credentials)
from bcbio.google.bc_metrics import (_write_to_worksheet,get_spreadsheet)

def write_run_report_to_gdocs(fc,qc,ssheet_title,encoded_credentials,wsheet_title=None):
    
    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title,encoded_credentials)
    if not client or not ssheet or qc is None:
        return False

    qc_metrics = qc.metrics()
    qc_stats = qc.getQCstats()
    run_cfg = qc.configuration()
    indexread = run_cfg.indexread()
    
    # Create the header row
    header = ["Lane","Description"]
    # Add the metric labels
    metric_lbl = []
    for metric in qc_metrics:
        if metric[0].endswith('_sd'):
            continue
        if metric[0] not in qc_stats:
            continue
        metric_lbl.append(metric[0])
        header.append(metric[1])
        
    # Iterate over the lanes of the flowcell and collect the data
    rows = []
    for lane in fc.get_lanes():
        # First the meta data
        row = [lane.get_name(),lane.get_description()]
        # Then the QC data
        for metric in metric_lbl:
            value = qc_stats[metric]
            sd_key = "%s_sd" % metric
            if sd_key in qc_stats:
                sd = qc_stats[sd_key]
            else:
                sd = None
            cell = ""
            for read in value.keys():
                if read.lstrip('read') in indexread:
                    continue
                val = "%s" % value[read][lane.get_name()]
                if sd is not None:
                    val += " +/- %s" % sd[read][lane.get_name()]
                if metric.find('cluster') >= 0:
                    cell = val
                    break
                cell += "%s: %s\n" % (read,val)
            row.append(cell)
        rows.append(row)
    
    if wsheet_title is None:
        wsheet_title = "%s_%s_%s" % (fc.get_fc_date(),fc.get_fc_name(),"QC")
    return _write_to_worksheet(client,ssheet,wsheet_title,rows,header,False)

                
        
    
    