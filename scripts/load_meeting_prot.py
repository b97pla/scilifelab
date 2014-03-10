"""Script to load info from trello and merge with info from latest meeting protocol."""

import bcbio.google
import bcbio.google.spreadsheet
import sys
import os
import operator
from optparse import OptionParser
from datetime import date, timedelta

CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')
credentials = bcbio.google.get_credentials({'gdocs_upload': {'gdocs_credentials': CREDENTIALS_FILE}})
CLIENT = bcbio.google.spreadsheet.get_client(credentials)

def get_ws(wsheet_title,ssheet):
    wsheet = bcbio.google.spreadsheet.get_worksheet(CLIENT,ssheet,wsheet_title)
    assert wsheet is not None, "Could not find worksheet %s within spreadsheet %s" % (wsheet_title,ssheet_title)
    content = bcbio.google.spreadsheet.get_cell_content(CLIENT,ssheet,wsheet)
    ws_key = bcbio.google.spreadsheet.get_key(wsheet)
    return ws_key, content

def get_meeting_info(content):
    "Feching info from old metting protocol"
    col=0
    FC = None
    data = {}
    for j,row in enumerate(content):
        if len(row[col].strip()) > 1:
            if row[col][1]=='.' or row[col][2]=='.':
                proj = row[col]
                data[proj] = {'info': row[col+1], 'flowcells':{}}
        elif len(row[col+1].strip()) > 1:
            if row[col+1][0]=='1':
                FC = row[col+1]
                data[proj]['flowcells'][FC] = [row[col+2]]
        elif not row[col+2].strip()=='' and FC:
            data[proj]['flowcells'][FC].append(row[col+2])
    return data

def sort_by_name(namelist):
    "Sorts dict alphabeticly by project sure name"
    name_dict={}
    for proj in namelist:
        name_dict[proj]=proj.split('.')[1].strip()
    sorted_name_dict = sorted(name_dict.iteritems(), key=operator.itemgetter(1))
    return sorted_name_dict

def pars_file(file, ongoing_deliveries, content):
    """parses output from script "update_checklist.py", which is loading runinfo
    from the trello board. ongoing_deliverues is a dict of info from old meeting
    that will be merged with the new info feched from the trello."""
    coming_deliveries = {}
    f=open(file,'r')
    content=f.readlines()
    nlines=len(content)
    data = coming_deliveries
    i=0
    proj=None
    while i < nlines:
        row = content[i]
        if row.strip() == "Ongoing":
            data = ongoing_deliveries
        if not row.strip()=='':
            isproj=row.split()[0]+'XXX'
            if isproj[1]=='.' or isproj[2]=='.':
                proj = row.split()[0].strip()
                if not proj in data.keys():
                    data[proj]={'info':' '.join(row.split()[1:]), 'flowcells':{}}
            elif proj and row.split()[0][0]!='1':
                if data[proj].has_key('info'):
                    if row.strip() not in data[proj]['info']:
                        data[proj]['info'] = data[proj]['info']+' - '+row
            elif proj and row.split()[0][0]=='1':
                if not data[proj]['flowcells'].has_key(row.strip()):
                    data[proj]['flowcells'][row.strip()]=[]
                proj=None
            i=i+1
        else:
            i=i+1
    return coming_deliveries, ongoing_deliveries

def update(sorted_names, col, info, ss_key, ws_key):
    """Uppdates the new meeting protocol with old and new info."""
    i=2
    for name in sorted_names:
        name = name[0]
        CLIENT.UpdateCell(i, col, name , ss_key, ws_key)
        CLIENT.UpdateCell(i, col+1, info[name]['info'] , ss_key, ws_key)
        fcs = sorted(info[name]['flowcells'].keys())
        for fc in fcs:
            i=i+1
            j=i
            CLIENT.UpdateCell(i, col+1, fc, ss_key, ws_key)
            for comment in info[name]['flowcells'][fc]:
                CLIENT.UpdateCell(j, col+2, comment, ss_key, ws_key)
                j=j+1
            if i!=j:
                i=j-1
        i=i+2

def main(old_wsheet, new_wsheet, file, ill_fix_it_myself):
    if os.path.exists("tempo.txt"): os.system("rm tempo.txt")
    if ill_fix_it_myself:
        old_wsheet = date.today()-timedelta(days=7)
        old_wsheet = old_wsheet.isoformat()
        new_wsheet = date.today().isoformat()
        file = "tempo.txt"
        os.system("update_checklist.py hugin-conf.yaml >> tempo.txt")
    ssheet_title='Genomics_Platform_Bioinformatics_Meeting_Deliveries'
    ssheet = bcbio.google.spreadsheet.get_spreadsheet(CLIENT,ssheet_title)
    assert ssheet is not None, "Could not find spreadsheet %s" % ssheet_title
    ss_key = bcbio.google.spreadsheet.get_key(ssheet)
    ws_key, content = get_ws(old_wsheet,ssheet)
    ongoing_deliveries = get_meeting_info(content)
    coming_deliveries, ongoing_deliveries = pars_file(file, ongoing_deliveries, content)
    sorted_ongoing_names = sort_by_name(ongoing_deliveries.keys())
    sorted_comming_names = sort_by_name(coming_deliveries.keys())
    ws_key, content = get_ws(new_wsheet,ssheet)
    update(sorted_ongoing_names, 1, ongoing_deliveries, ss_key, ws_key)
    update(sorted_comming_names, 4, coming_deliveries,ss_key, ws_key)
    if os.path.exists("tempo.txt"): os.system("rm tempo.txt")

if __name__ == "__main__":
    parser = OptionParser(usage = "analysis_reports.py <project id> [Options]")
    parser.add_option("-o", "--old_wsheet", dest="old_wsheet",default=None,
    help = "Get old notes from this wsheet")
    parser.add_option("-n", "--new_wsheet", dest="new_wsheet",default=None,
    help = "Load this wsheet with old notes and new notes.")
    parser.add_option("-f", "--file_dump", dest="file",default=None,
    help = "Get info about new runs from this file. Give the whole filepath")
    parser.add_option("-i", "--ill_fix_it_myself", dest="ill_fix_it_myself",action="store_true", default=None,
    help = """Will try to run update_checklist.py, will search for old wsheet named as the
    iso date one week from now, and will load a new wsheet named by todays iso date if existing""")
    (options, args) = parser.parse_args()

    main(options.old_wsheet, options.new_wsheet, options.file, options.ill_fix_it_myself)


