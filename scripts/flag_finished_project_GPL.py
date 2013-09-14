#!usr/bin/env -python

import re
import os
import base64
import argparse
import datetime
import subprocess
import gdata
import gdata.docs.service
import gdata.spreadsheet.service
from scilifelab.google import get_credentials
from scilifelab.google.google_docs import SpreadSheet

## declaring arguments and values for the script
parser = argparse.ArgumentParser()
parser.add_argument("Config_file", type=str, help="Path to the config file that has info for the run")
args = parser.parse_args()
fCnt = 0
search1 = 'all raw-data delivered.*sign'
search2 = 'all raw-data delivered.*date'
search3 = 'project name'
config = {}
t = str(datetime.datetime.utcnow().isoformat()).split('.')[0]

##looks for config file and gets the values needed for the run
config_file = args.Config_file
if os.path.exists(config_file):
	file_config = open(config_file,'r')
	for line in file_config.readlines():
		vals = line.strip().split(':')
		config[vals[0]] = vals[1]
else:
	print "\nConfig file doesn't exist, Provide valid file\n"
	raise SystemExit

print "\nProgram started at {}".format(t) ##log
print "--------------------------------------\n" ##log
print "Fetching Spreadsheet from google drive..." ##log

## method to merge two header rows to one
def get_head(head1,head2):
	head_main = []
	hd1 = [h.replace('\n','') for h in head1]
	hd2 = [h.replace('\n','') for h in head2]
	for i,h in enumerate(hd1):
		if h != '':
			head_main.append("{}_{}".format(h,hd2[i]))
			pre = h
		elif h == '':
			head_main.append("{}_{}".format(pre,hd2[i]))
	return head_main;

## method the check if already the project is been flagged
def check_flag(p_path):
	sam_cnt, sam_nm = (0, [])
	## gets only the directories i.e. samples from project folder ##
	sample_dir = [h for h in os.listdir(p_path) if not os.path.isfile("{}/{}".format(p_path,h))]
	for sample in sample_dir:
		flag = os.path.join(p_path,sample,'FINISHED_AND_DELIVERED')
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

## getting the spreadsheet from given name and appropriate worksheet
credentials = get_credentials(config['gCred_file'])
ssheet = SpreadSheet(credentials,config['ssheet_nm'])
wksheet = ssheet.get_worksheet(config['wksheet_nm'])
cell_feed = ssheet.get_cell_feed(wksheet)

print "**** DONE ****\n" ##log
print "Parsing worksheet and obtaining information..." ##log

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
print "{} projects are marked as deliverd in GPL".format(str(len(projects_done))) ##log
print "Starting to touch-finished for unflagged projects\n" ##log
print "Project\tNum_of_Samples\tDelivered_date\tDelivered_by" ##log

## from the obtained list project unflagged are flagged
upp_path = config['upp_path']
for p in projects_done:
	proj, person, deliv_date = p
	uP = os.path.join(upp_path,proj)
	if os.path.exists(uP):
		flag_status,samples_unmarked = check_flag(uP)
		if flag_status != 'full_marked':
			for s in samples_unmarked:
				cmd = subprocess.check_call(['pm','production','touch-finished',proj,'-S',s,'--force','--quiet'])
			fCnt = fCnt+1
			print "{}\t{}\t{}\t{}".format(proj,str(len(samples_unmarked)),deliv_date,person) ##log
	else:
		print "RUN TERMINATIN: project directory not found for {}".format(proj) ##log
		raise SystemExit

## log
print "\n*** Touch done ****"
print "\nRUN SUMMARY:"
print "------------"
print "Project marked delivered in GPL: {}".format(str(len(projects_done)))
print "No. of project already flagged: {}".format(str(len(projects_done)-fCnt))
print "No. of project flagged in this run: {}\n".format(str(fCnt))
