#!/usr/bin/env python
"""Functions for getting barcode statistics from demultiplexing.
"""

import copy
import bcbio.google.connection
import bcbio.google.document
from bcbio.google import spreadsheet as g_spreadsheet
from bcbio.google import _to_unicode
from bcbio.log import logger2
import bcbio.solexa.flowcell
import bcbio.pipeline.flowcell

# The structure of the demultiplex result
BARCODE_STATS_HEADER = [
         ['Project name', 'project_name'],
         ['Lane', 'lane'],
         ['Lane description', 'description'],
         ['Sample name', 'sample_name'],
         ['bcbb internal barcode index', 'bcbb_barcode_id'],
         ['Barcode name', 'barcode_name'],
         ['Barcode sequence', 'barcode_sequence'],
         ['Barcode type', 'barcode_type'],
         ['Demultiplexed read (pair) count', 'read_count'],
         ['Demultiplexed read (pair) count (millions)', 'rounded_read_count'],
         ['Comment', 'comment']
        ]

# The structure of the sequencing result
SEQUENCING_RESULT_HEADER = [
                 ['Sample name', 'sample_name'],
                 ['Run', 'run'],
                 ['Lane', 'lane'],
                 ['Read (pair) count', 'read_count'],
                 ['Read (pair) count (millions)', 'rounded_read_count'],
                 ['Barcode sequence', 'barcode_sequence'],
                 ['Comment', 'comment'],
                 ['Pass', 'pass']
                ]


def _create_header(header, columns):
    header = copy.deepcopy(header)
    names = []
    for column in columns:
        for i, head in enumerate(header):
            if (len(head) > 1 and head[1] == column):
                names.append(head[0])
                header.pop(i)
    for head in header:
        names.append(head[0])
    return names

def _header_index(header, label, header2=[]):
    """Returns the index of the label in the supplied header (assumed to have a structure of [[header, label]]).
       or -1 if not found. If a second header is supplied, will instead return the index of the corresponding 
       header string in that.
    """
    try:
        ind = [head[1] for head in header].index(label)
    except:
        return -1
        
    if len(header2) > 0:
        try:
            ind = header2.index(header[ind][0])
        except:
            return -1
    
    return ind
        
    
def get_spreadsheet(ssheet_title, encoded_credentials):
    """Connect to Google docs and get a spreadsheet"""

    # Convert the spreadsheet title to unicode
    ssheet_title = _to_unicode(ssheet_title)

    # Create a client class which will make HTTP requests with Google Docs server.
    client = g_spreadsheet.get_client()
    bcbio.google.connection.authenticate(client, encoded_credentials)

    # Locate the spreadsheet
    ssheet = g_spreadsheet.get_spreadsheet(client, ssheet_title)

    # Check that we got a result back
    if not ssheet:
        logger2.warn("No document with specified title '%s' found in \
                      GoogleDocs repository" % ssheet_title)
        return (None, None)

    return (client, ssheet)


def _write_project_report_to_gdocs(client, ssheet, flowcell):
    """Writes report data to corresponding worksheets in a google docs
    spreadsheet.
    """
    # Flatten the project_data structure into a list
    samples = {}
    for sample in flowcell.get_samples():
        if sample.get_name() in samples:
            samples[sample.get_name()].add_sample(sample)
        else:
            samples[sample.get_name()] = sample

    column_headers = [col_header[0] for col_header in SEQUENCING_RESULT_HEADER[:6]]

    success = True
    for sample in samples.values():
        run_name = "{}_{}".format(flowcell.fc_date, flowcell.fc_name)
        wsheet_title = run_name
        row = (sample.get_name(), \
               run_name, \
               sample.get_lane(), \
               sample.get_read_count(), \
               sample.get_rounded_read_count(), \
               sample.barcode_sequence)

        success &= _write_to_worksheet(client, \
                                       ssheet, \
                                       wsheet_title, \
                                       [row], \
                                       column_headers, \
                                       append=True,
                                       keys=[head[0] for head in SEQUENCING_RESULT_HEADER if head[1] in ["sample_name", "lane", "barcode_sequence"]])

    return success


def write_project_report_summary_to_gdocs(client, ssheet):
    """Summarize the data from the worksheets and write them to a "Summary"
    worksheet.
    """
    # Summary data
    flowcells = {}
    samples = {}
    # Get the list of worksheets in the spreadsheet
    wsheet_feed = g_spreadsheet.get_worksheets_feed(client, ssheet)
    # Loop over the worksheets and parse the data from the ones that contain
    # flowcell data.

    for wsheet in wsheet_feed.entry:
        wsheet_title = wsheet.title.text
        if wsheet_title.endswith("_QC"):
            continue
        
        # Use the bcbio.solexa.flowcell.get_flowcell_info method to determine if the wsheet title contains a valid flowcell id
        try:
            bcbio.solexa.flowcell.get_flowcell_info(wsheet_title)
        except ValueError:
            continue

        # Get the worksheet header
        wsheet_header = g_spreadsheet.get_header(client, ssheet, wsheet)

        wsheet_data = g_spreadsheet.get_cell_content(client, ssheet, wsheet, '2')
        delim = ';'
        # Map the column names to the correct index using the header
        sample_col, run_col, lane_col, count_col, bc_col = [_header_index(SEQUENCING_RESULT_HEADER,col_name,wsheet_header) for col_name in ('sample_name', 'run', 'lane', 'read_count', 'barcode_sequence')]
        
        # Add the results from the worksheet to the summarized data.
        for row in wsheet_data:

            sample_name, run_name, lane_name, read_count, barcode_sequence = [row[col] if col >= 0 else None for col in (sample_col, run_col, lane_col, count_col, bc_col)]
                 
            data = {"name": sample_name,
                    "read_count": read_count,
                    "sequence": barcode_sequence}

            lane = bcbio.pipeline.flowcell.Lane({"lane": lane_name})
            sample = bcbio.pipeline.flowcell.BarcodedSample(data, lane)

            if sample_name in samples:
                samples[sample_name]["object"].add_sample(sample, delim)
                samples[sample_name]["flowcells"] += "{}{}".format(delim, run_name)
                if not samples[sample_name]["object"].barcode_sequence and barcode_sequence:
                    samples[sample_name]["object"].barcode_sequence = barcode_sequence

            else:
                samples[sample_name] = {"object": sample, "flowcells": run_name}

    wsheet_title = "Summary"

    # Try getting already existing 'comment' and 'pass' values.
    name_data = {}
    existing_summary_wsheet = g_spreadsheet.get_worksheet(client, ssheet, wsheet_title)
    if existing_summary_wsheet:
        summary_header = g_spreadsheet.get_header(client, ssheet, existing_summary_wsheet)
        summary_data = g_spreadsheet.get_cell_content(client, ssheet, existing_summary_wsheet, '2')
        sample_col, comment_col, pass_col = [_header_index(SEQUENCING_RESULT_HEADER,col_name,summary_header) for col_name in ('sample_name', 'comment', 'pass')]
        
        for content in summary_data:
            sample_name, comment, pass_field = [content[col] if col >= 0 else None for col in (sample_col, comment_col, pass_col)]
            name_data[sample_name] = [comment, pass_field]

    # Flatten the project_data structure into a list
    rows = []
    for sample_data in samples.values():
        sample = sample_data["object"]
        flowcells = sample_data["flowcells"]

        sample_name = sample.get_name()
        comment, pass_field = name_data.get(sample_name, ["", ""])

        row = [sample_name, \
               flowcells, \
               sample.get_lane(), \
               sample.get_read_count(), \
               sample.get_rounded_read_count(), \
               sample.barcode_sequence, \
               comment, \
               pass_field]

        rows.append(row)

    # Write the data to the worksheet
    column_headers = [col_header[0] for col_header in SEQUENCING_RESULT_HEADER]
    return _write_to_worksheet(client, \
                               ssheet, \
                               wsheet_title, \
                               rows, \
                               column_headers, \
                               False)


def write_run_report_to_gdocs(fc, fc_date, fc_name, ssheet_title, \
    encoded_credentials, wsheet_title=None, append=False, split_project=False):
    """Upload the barcode read distribution for a run to google docs.
    """

    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title, encoded_credentials)
    if not client or not ssheet:
        return False

    # Get the projects in the run
    projects = fc.get_project_names()
    logger2.info("Will write data from the run {}_{} for " \
        "projects: {!r}".format(fc_date, fc_name, "', '".join(projects)))

    # If we will split the worksheet by project, use the project
    # names as worksheet titles.
    success = True
    header = _create_header(BARCODE_STATS_HEADER, fc.columns())
    if split_project:
        # Filter away the irrelevent project entries and write the
        # remaining to the appropriate worksheet.
        for project in projects:
            pruned_fc = fc.prune_to_project(project)
            success &= _write_to_worksheet(client, ssheet, project, \
                pruned_fc.to_rows(), header, append, \
                keys=[head[0] for head in BARCODE_STATS_HEADER if head[1] in ["sample_name", "lane", "barcode_sequence"]])

    # Otherwise, set the default title of the worksheet to be a string of
    # concatenated date and flowcell id.
    else:
        if wsheet_title is None:
            wsheet_title = "{}_{}".format(fc_date, fc_name)

        success &= _write_to_worksheet(client, ssheet, wsheet_title, \
            fc.to_rows(), header, append, \
            keys=[head[0] for head in BARCODE_STATS_HEADER if head[1] in ["sample_name", "lane", "barcode_sequence"]])

    return success


def _write_to_worksheet(client, ssheet, wsheet_title, rows, header, append, keys=[]):
    """Generic method to write a set of rows to a worksheet on google docs.
    """
    # Convert the worksheet title to unicode
    wsheet_title = _to_unicode(wsheet_title)

    # Add a new worksheet, possibly appending or replacing a pre-existing
    # worksheet according to the append-flag.
    wsheet = g_spreadsheet.add_worksheet(client, \
                                         ssheet, \
                                         wsheet_title, \
                                         len(rows) + 1, \
                                         len(header), \
                                         append)
    if wsheet is None:
        logger2.error("ERROR: Could not add a worksheet {!r} to " \
            "spreadsheet {!r}".format(wsheet_title, ssheet.title.text))
        return False
    
    # If keys are specified (will correspond to indexes in the header), delete pre-existing rows with matching keys
    if append and len(keys) > 0:
        wsheet_data = g_spreadsheet.get_cell_content(client, ssheet, wsheet, '2')
        wsheet_header = g_spreadsheet.get_header(client, ssheet, wsheet)
        try:
            wsheet_indexes = [wsheet_header.index(key) for key in keys]
            header_indexes = [header.index(key) for key in keys]
        except ValueError:
            logger2.warn("WARNING: Could not identify correct header for duplicate detection")
        else:
            for row in rows:
                try:
                    key = "#".join([row[i] for i in header_indexes])        
                    for i, wrow in enumerate(wsheet_data):
                        wkey = "#".join([wrow[j] for j in wsheet_indexes])
                        if wkey == key:
                            g_spreadsheet.delete_row(client, ssheet, wsheet, i+1)
                            wsheet_data.pop(i)
                            break
                except:
                    logger2.warn("WARNING: Could not identify/replace duplicate rows")

    # Write the data to the worksheet
    success = g_spreadsheet.write_rows(client, ssheet, wsheet, header, rows)
    if success:
        logger2.info("Wrote data to the {!r}:{!r} " \
                     "worksheet".format(ssheet.title.text, wsheet_title))
    else:
        logger2.error("ERROR: Could not write data to the {!r}:{!r} " \
                      "worksheet".format(ssheet.title.text, wsheet_title))
    return success
