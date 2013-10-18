
import os
import datetime
import re
import string
from cStringIO import StringIO
from scilifelab.bcbio.qc import FlowcellRunMetricsParser
import scilifelab.log
import scilifelab.google as google
from scilifelab.google.google_docs import SpreadSheet

LOG = scilifelab.log.minimal_logger(__name__)

def _column_mapping():
    """Return the mapping between casava report headers and google docs headers"""
    mapping = [["Project", "Project name"],
               ["Sample ID","Sample name"],
               ["Description", "Description"],
               ["Lane", "Lane"],
               ["", "Read pair count"],
               ["", "Read pairs (Mbases)"],
               ["# Reads", "Total reads"],
               ["Index", "Barcode sequence"],
               ["% of >= Q30 Bases (PF)", "% of >= Q30 Bases (PF)"],
               ["Mean Quality Score (PF)", "Mean Quality Score (PF)"],
               ["% Perfect Index Reads", "% Perfect Index Reads"],
               ["% One Mismatch Reads (Index)", "% One Mismatch Reads (Index)"],
               ["Yield (Mbases)", "Yield (Mbases)"],
               ["Sample Ref", "Genome build"],
               ["", "Run setup"],
               ["% PF", "% Pass filter"],
               ["% of raw clusters per lane", "% of raw clusters per lane"],
               ["Control", "Control sample (Y/N)"]]
    return mapping

def _column_header():
    return [m[1] for m in _column_mapping() if len(m[1]) > 0]

def _demultiplex_spreadsheet(run_date):
    """Get the appropriate demultiplex spreadsheet according to the date of the run
    """
    ts = datetime.datetime.strptime(run_date,"%y%m%d")
    quarter = 1+(ts.month-1)/3
    return "Demultiplex Counts {} Q{}".format(str(ts.year),str(quarter))

def _run_setup(reads):
    read_cycles = [r.get('NumCycles','?') for r in reads if r.get('IsIndexedRead','N') == 'N']
    index_cycles = [r.get('NumCycles','?') for r in reads if r.get('IsIndexedRead','N') == 'Y']
    nreads = len(read_cycles)
    nindex = len(index_cycles)
    
    c = list(set(read_cycles))
    if len(c) == 1:
        setup = "{}x{}".format(str(nreads),c[0])
    else:
        setup = ",".join(c)
    
    return setup

def upload_to_gdocs(fcdir, credentials_file=None, gdocs_folder=None):

    output_data = {'stdout':StringIO(), 'stderr':StringIO(), 'debug':StringIO()}
    
    if not os.path.exists(fcdir):
        LOG.error("The run folder, {} does not exist!".format(os.path.basename(fcdir)))
        return output_data
    
    credentials = google.get_credentials(credentials_file)
    if credentials is None:
        LOG.error("Could not parse the Google Docs credentials")
        return output_data
    
    metrics = collect_metrics(fcdir)
    samples = _format_samples(metrics)
    
    ssheet_name = _demultiplex_spreadsheet(metrics['RunInfo'].get('Date',None))
    ssheet = SpreadSheet(credentials,ssheet_name)
    ssheet.move_to_folder(gdocs_folder)
    
    run_id = metrics['RunInfo']['Id'].split("_")
    wsheet_name = "_".join([run_id[0],run_id[-1]])
    
    # Write the metrics for the entire flowcell
    write_flowcell_metrics(samples, ssheet, wsheet_name)
    
    # Write project-centered metrics
    projects = list(set([sample.get('Project name','') for sample in samples]))
    for project in projects:
        if project in ['Undetermined_indices','']:
            continue
        project_samples = [sample for sample in samples if sample.get('Project name','') == project]
        # Insert the run name as description
        for sample in project_samples:
            sample['Description'] = wsheet_name
            
        ssheet_name = "{}_sequencing_results".format(project)
        ssheet = SpreadSheet(credentials,ssheet_name)
        ssheet.move_to_folder(gdocs_folder)
        # Truncate the summary worksheet so that it won't show the wrong information in case upload fails
        write_flowcell_metrics([], ssheet, "Summary")
        project_samples = summarize_project(ssheet,{wsheet_name: project_samples})
        write_flowcell_metrics(project_samples, ssheet, wsheet_name)
        
        # Create the summary over all worksheets in the project
        summary_samples = summarize_project(ssheet)
        write_flowcell_metrics(summary_samples, ssheet, "Summary")
    
    return output_data

def summarize_project(ssheet, project_data=None):
    """Parse the flowcell worksheets for a project and sum the counts for each sample
    """
    
    # If the data to be summarized was not supplied, get the feed of all worksheets
    if project_data is None:
        project_data = {}
        for wsheet in ssheet.get_worksheets_feed().entry:
            wsheet_title = wsheet.title.text
            
            # Skip a QC worksheet
            if wsheet_title.endswith("_QC"):
                continue
            
            # Skip a worksheet if it doesn't have the [DATE]_[FCID] name pattern
            s = wsheet_title.split("_")
            if len(s) != 2:
                continue
            try:
                datetime.datetime.strptime(s[0],'%y%m%d')
            except ValueError:
                continue
            
            # Get the worksheet contents
            wsheet_data = ssheet.get_cell_content(wsheet)
            project_data[wsheet_title] = []
            
            # Convert the data into a dictionary
            for row in wsheet_data[1:]:
                project_data[wsheet_title].append(dict(zip(wsheet_data[0],row)))
        
    summary = {}
    header = _column_header()
    
    for wsheet in sorted(project_data.keys()):
        for n, sample in enumerate(project_data[wsheet]):
            sample_key = _match_key_name(summary.keys(), sample['Sample name'])
            sample_summary = summary.get(sample_key,{})
            
            if len(sample_key) > sample['Sample name']:
                sample_summary['Sample name'] = sample_key
            else:
                sample_summary['Sample name'] = sample['Sample name']
                
            if 'Project name' in sample:
                sample_summary['Project name'] = sample['Project name']
            
            sample_summary['Description'] = sample_summary.get('Description',"").split(";") + [sample.get('Description','')]
            sample_summary['Description'] = ";".join([s for s in sample_summary['Description'] if len(s) > 0])
            
            sample_summary['Lane'] = sample_summary.get('Lane',"").split(";") + [sample.get('Lane','')]
            sample_summary['Lane'] = ";".join([s for s in sample_summary['Lane'] if len(s) > 0])
            
            sample_summary['Total reads'] = sample_summary.get('Total reads',0) + int(sample.get('Total reads',sample.get('Read (pair) count',0)))
            sample_summary['Read pair count'] = sample_summary.get('Read pair count',0) + int(sample.get('Read pair count',sample.get('Read (pair) count',0)))
            sample_summary['Read pairs (Mbases)'] = sample_summary.get('Read pairs (Mbases)',0) + float(sample.get('Read pairs (Mbases)',sample.get('Read (pair) count (millions)',0)))
            
            sample_summary['Barcode sequence'] = sample_summary.get('Barcode sequence',"").split(";") + [sample.get('Barcode sequence','')]
            sample_summary['Barcode sequence'] = ";".join([s for s in sample_summary['Barcode sequence'] if len(s) > 0])
    
            summary[sample_key] = sample_summary
            
    return summary.values()

            
def _match_key_name(haystack, needle):
    
    min_match = 5
    for straw in haystack:
        if straw == needle:
            return straw
        if straw.startswith(needle) and len(needle) >= min_match:
            return straw
        if needle.startswith(straw) and len(straw) >= min_match:
            return straw
    return needle
        
def _format_samples(metrics):
    
    # Get the casava - gdocs mapping and the header columns
    mapping = _column_mapping()
    
    # Determine the run setup
    reads = metrics['RunInfo'].get('Reads')
    run_setup = _run_setup(reads)
    
    samples = []
    for sample in metrics.get('Barcode_lane_statistics',[]):
        sample_data = {c[1]: sample.get(c[0],"") for c in mapping if len(c[0]) > 0}
        
        # Add the run setup
        sample_data['Run setup'] = run_setup
        
        # Calculate the number of read pairs
        total = int(sample_data['Total reads'].replace(",",""))
        pairs = total/2
        if len([r for r in reads if r.get('IsIndexedRead','N') == 'N']) == 1:
            pairs = total
        sample_data['Total reads'] = str(total)
        sample_data['Read pair count'] = str(pairs)
        sample_data['Read pairs (Mbases)'] = str(round(pairs/1000000.,2))
        
        # Format some of the columns a bit nicer
        sample_data['Project name'] = sample_data['Project name'].replace("__",".")
        sample_data['Description'] = sample_data['Description'].replace("__",".") 
        sample_data['Yield (Mbases)'] = sample_data['Yield (Mbases)'].replace(",","")
        samples.append(sample_data)
        
    # Append the undemultiplexed index numbers
    for lane in sorted(metrics.get('Undemultiplexed',{}).keys()):
        data = metrics['Undemultiplexed'][lane].get('undemultiplexed_barcodes',{})

        for i in range(len(data.get('count',[]))):
            sample_data = {'Description': data['index_name'][i],
                           'Barcode sequence': data['sequence'][i],
                           'Total reads': data['count'][i],
                           'Lane': lane,
                           'Sample name': 'Undemultiplexed index'}
            
            # Calculate the number of read pairs
            total = int(sample_data['Total reads'])
            pairs = total/2
            if len([r for r in reads if r.get('IsIndexedRead','N') == 'N']) == 1:
                pairs = total
            sample_data['Total reads'] = str(total)
            sample_data['Read pairs (Mbases)'] = str(round(pairs/1000000.,2))
            samples.append(sample_data)
    
    return samples

def write_flowcell_metrics(samples, ssheet, wsheet_name):
    """Write the supplied metrics for the flowcell to Google Docs"""
    
    # Get the casava - gdocs mapping and the header columns
    header = _column_header()
    
    # Populate the rows
    rows = []
    for sample in samples:
        rows.append([sample.get(h,"") for h in header])
        
    # Create the target worksheet
    wsheet = ssheet.add_worksheet(wsheet_name, rows=len(rows)+1, cols=len(header))
    return ssheet.write_rows(wsheet,header,rows)


def collect_metrics(path):
    parser = FlowcellRunMetricsParser(path)
    run_info = parser.parseRunInfo()
    fcid = run_info.get('Flowcell',None)
    if fcid is None:
        LOG.error("Could not parse flowcell id from RunInfo.xml")
        return {}
    
    # Insert a dummy character as the parse method expects a flowcell position
    metrics = parser.parse_demultiplex_stats_htm(fcid)
    metrics['RunInfo'] = run_info
    
    # Get the undemultiplexed indexes
    undemux = parser.parse_undemultiplexed_barcode_metrics(fcid)
    metrics['Undemultiplexed'] = undemux
    
    return metrics
