#!usr/bin/env -python

import gdata
import gdata.docs.service
import gdata.spreadsheet.service
import base64
import re
import os
from os.path import expanduser
import datetime
import subprocess
import sys

##script called time
t = str(datetime.datetime.now()).split('.')[0]

print "\nProgram started at "+t ##log
print "--------------------------------------\n" ##log

##looking for gdoc credentials file in opt/config
gdoc = sys.argv[1]
try:
	gdoc_cred = open(gdoc,'r').readline().strip()
	gdoc_mail, gdoc_pass = base64.b64decode(gdoc_cred).split(':')
except IOError:
	print "\nRUN TERMINATING: No gdoc credential file found in \"opt/config/\"\n"
	

## method to merge two header rows to one
def get_head(head1,head2):
	head_main = []
	hd1 = [h.replace('\n','') for h in head1]
	hd2 = [h.replace('\n','') for h in head2]
	for i,h in enumerate(hd1):
		if h != '':
			head_main.append(h+'_'+hd2[i])
			pre = h
		elif h == '':
			head_main.append(pre+'_'+hd2[i])
	return head_main;

## method the check if already the project is been flagged
def check_flag(p_path):
	sam_cnt, sam_nm = (0, [])
	## gets only the directories i.e. samples from project folder ##
	sample_dir = [h for h in os.listdir(p_path) if not os.path.isfile(p_path+'/'+h)]
	for sample in sample_dir:
		flag = p_path+'/'+sample+'/'+'FINISHED_AND_DELIVERED'
		if not os.path.exists(flag):
			sam_cnt = sam_cnt+1
			sam_nm.append(sample)
	if sam_cnt == 0:
		status = 'full_marked'
	elif sam_cnt == len(sample_dir):
		status = 'none_marked'
	else:
		status = 'some_marked'
	return status,sam_nm;

## setting up client document for google docs
client = gdata.spreadsheet.service.SpreadsheetsService()
client.email = gdoc_mail
client.password = gdoc_pass
client.source = "Touching finished projects"
client.ProgrammaticLogin()

## setting key and title for the spreadsheet and worksheet
sp_key = sys.argv[2]
wk_nm = 'Ongoing'
search1 = 'all raw-data delivered.*sign'
search2 = 'all raw-data delivered.*date'
search3 = 'project name'
fCnt = 0

print "Fetching worksheet from google drive..." ##log

## making a query and retrieving the worksheet with exact name
wP = {'title':wk_nm, 'title-exact':'True'}
wQ = gdata.spreadsheet.service.DocumentQuery(params=wP)
wk_sheets = client.GetWorksheetsFeed(key=sp_key, query=wQ)

## second check to make sure proceeding with correct worksheet
for wk in wk_sheets.entry:
	if wk.title.text == wk_nm:
		wk_main = wk
		break

print "**** DONE ****\n" ##log
print "Parsing the cell feed and collecting information..." ##log

## setting cell query and getting cell's feed as object
wk_id = wk_main.id.text.split('/')[-1]
rw_st, rw_en = ('1', wk_main.row_count.text)
col_st, col_en =  ('1', wk_main.col_count.text)
cell_P = {'min-row': rw_st, 'max-row': rw_en,
          'min-col': col_st, 'max-col': col_en,
          'return-empty': 'True'
          }
cell_Q = gdata.spreadsheet.service.CellQuery(params=cell_P)
cell_feed = client.GetCellsFeed(key=sp_key, wksht_id=wk_id, query=cell_Q)

## iterating through cell's content and gives list of prject signed delivered
col_tot = int(cell_feed.col_count.text)
projects_done, cell_content, row, ck, rNum = ([], [], [], 0, 1)
for cell in cell_feed.entry:
	row.append(cell.content.text or "")
	ck = ck+1
	if ck == col_tot:
		if rNum == 1:
			head1 = row
		elif rNum == 2:
			head2 = row
			header = get_head(head1,head2)
			status_sign = [nm for nm in header if re.search(search1,nm,re.I)][0]
			status_date = [nm for nm in header if re.search(search2,nm,re.I)][0]
			proj_col = [nm for nm in header if re.search(search3,nm,re.I)][0]
			sign_ind = header.index(status_sign)
			date_ind = header.index(status_date)
			proj_ind = header.index(proj_col)
		else:
			if row[sign_ind] != '':
				projects_done.append([row[proj_ind],row[sign_ind],row[date_ind]])
		rNum = rNum+1
		row, ck = ([], 0)

print "**** DONE ****\n" ##log	
print str(len(projects_done))+" projects are marked as deliverd in GPL" ##log
print "Starting to touch-finished for unflagged projects\n" ##log
print "Project\tNum_of_Samples" ##log

## from the obtained list project unflagged are flagged
fold_path = 'path/for/uppmax/projects/here/'
for p in projects_done:
	pd = fold_path+'/'+p
	if os.path.exists(pd):
		flag_status,samples_unmarked = check_flag(pd)
		if flag_status != 'full_marked':
			for s in samples_unmarked:
				subprocess.call('pm production touch-finished '+p+' -S '+s+' --force --quiet', shell=True)
			fCnt = fCnt+1
			print p+'\t'+str(len(samples_unmarked)) ##log
	else:
		print "RUN TERMINATIN: project directory not found for "+p

##log
print "\n*** Touch done ****"
print "\nRUN SUMMARY:"
print "------------"
print "Project marked delivered in GPL: "+str(len(proj_ls))
print "No. of project already flagged: "+str(len(proj_ls)-fCnt)
print "No. of project flagged in this run: "+str(fCnt)+"\n"
