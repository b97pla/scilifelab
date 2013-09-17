import bcbio.google
import bcbio.google.spreadsheet
import sys
import os

if len(sys.argv) < 3:
        print """
This script loads Reads after duplicate removal into the 20158 table given that 
you have run the analysis pipeline so that you have a stat.json file with mapping 
statistics in the analysis directory.
It will overwrite whatever is writen in the column "Reads after duplicate removal" 
in table 20158. You need to check what VERSION of the document that you what to upload to.

IMPORTANT: Version 2 of the document 20158 has a mistake in the column named: 
           "Total number of reads after duplicate removal (Millions)".
           This mistake should vave been fixed for most projects, but if you see
           that the fielsd of this column states "VALUE", please contact Maya 


Usage:
To be run from the analysis directory, eg: 

/proj/a2012043/private/nobackup/projects/t_olsson_11_02/intermediate/20120503A_hiseq2000

load_dupRem_20158.py <project id> <version> -s

        <project id>    eg: T.Olsson_11_02
        <version>       document version of 20158 integer - either 1,2 or 3 
        [-s]            optional for single end 
        """
        sys.exit()

PROJECT_ID      = sys.argv[1]
VERSION         = sys.argv[2]
try:
        if sys.argv[3] ==  '-s':
                single = True
        else:
                single=False
except:
        single=False
        pass

CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')

dup_rem_col = "Total number of reads after duplicate removal (Millions)"
scilife_col = "Sample name (SciLifeLab)"


if VERSION == "2":
        ssheet_title    = "%s_20158_02 QC for HiSeq sequencing results" % (PROJECT_ID)
elif VERSION == "3":
        ssheet_title    = "%s_20158_03 QC for HiSeq sequencing results" % (PROJECT_ID)
elif VERSION == "5":
        ssheet_title    = "%s_20158_05 QC for HiSeq sequencing results" % (PROJECT_ID)
        dup_rem_col     = "Total number after duplicate removal"
        scilife_col = "Sample"
elif VERSION == "6":
        ssheet_title    = "%s_20158_06 QC for HiSeq sequencing results" % (PROJECT_ID)
        dup_rem_col     = "Total number after duplicate removal"
        scilife_col = "Sample"
else:
        sys.exit("Unknown VERSIONnumber")

wsheet_title = 'Sheet1'

try:
        dict = eval(open('stat.json').read())
except:
        sys.exit('could not fine mapping statistics file stat.json')

#-------------
credentials = bcbio.google.get_credentials({'gdocs_upload': {'gdocs_credentials': CREDENTIALS_FILE}})
client = bcbio.google.spreadsheet.get_client(credentials)

ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
assert ssheet is not None, "Could not find spreadsheet %s" % ssheet_title

wsheet = bcbio.google.spreadsheet.get_worksheet(client,ssheet,wsheet_title)
assert wsheet is not None, "Could not find worksheet %s within spreadsheet %s" % (wsheet_title,ssheet_title)

content = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet)

ss_key = bcbio.google.spreadsheet.get_key(ssheet)
ws_key = bcbio.google.spreadsheet.get_key(wsheet)
#--------------

reads_colindex = 0
names_colindex = 0
rowindex = 0
for j,row in enumerate(content):
        if (names_colindex == 0) | (reads_colindex == 0):
                for i,col in enumerate(row):
                        if str(col).strip() == scilife_col:
                                names_colindex = i+1
                        if str(col).strip() == dup_rem_col:
                                reads_colindex = i+1
        name = str(row[names_colindex-1]).strip()
        try:
                name=name.split('_')[0]+'_'+name.split('_')[1]
        except:
                pass
        if dict.has_key(name):
                if single:
                        print dict[name]['aft_dup_rem']['Uniquely mapped']
                        M_reads_aft_dup_rem=str(round(float(dict[name]['aft_dup_rem']['Uniquely mapped'])/1000000.0,2))
                else:
                        R1=dict[name]['aft_dup_rem']['Read-1']
                        R2=dict[name]['aft_dup_rem']['Read-2']
                        M_reads_aft_dup_rem=str(round((float(R2)+float(R1))/2000000,2))
                client.UpdateCell(j+1, reads_colindex, M_reads_aft_dup_rem, ss_key, ws_key)
                print name+' '+M_reads_aft_dup_rem











