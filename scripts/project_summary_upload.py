#!/usr/bin/env python
import sys
import os
import time
from  datetime  import  datetime
from uuid import uuid4
import hashlib
from optparse import OptionParser
import logging

import bcbio.google
import scilifelab.google.project_metadata as pmeta
import bcbio.pipeline.config_loader as cl
from bcbio.google import _to_unicode, spreadsheet
import couchdb

def get_proj_inf(WS_projects,project_name_swe, samp_db, proj_db, client, config):
	project_name = _replace_ascii(_to_unicode(project_name_swe))
	key = find_proj_from_view(proj_db, project_name)
	if not key: key = uuid4().hex

        logger.info('Handling proj %s %s' % (project_name, key))

        obj={'application': '',
	     'customer_reference': '',
	     'min_m_reads_per_sample_ordered': '',
	     'no_of_samples': '',
             'entity_type': 'project_summary',
             'uppnex_id': '',                 
             'samples': {},
             'project_name': project_name, 
	     'project_id':'',
             '_id': key}


	### Genomics Project list
	p = pmeta.ProjectMetaData(project_name, config)
	if p.project_name is None:
		p = pmeta.ProjectMetaData(project_name_swe, config)
	if p.project_name is None:
		logger.warning('Google Document Genomics Project list: %s not found' % project_name) 
	else:
		if p.min_reads_per_sample.strip() != '':
                	try: obj['min_m_reads_per_sample_ordered'] = round(float(p.min_reads_per_sample),2)
			except: obj['min_m_reads_per_sample_ordered'] = p.min_reads_per_sample
		if p.no_samples.strip() != '':
			try:obj['no_of_samples'] = int(p.no_samples)
			except:obj['no_of_samples'] = p.no_samples
                obj['uppnex_id'] = p.uppnex_id
		obj['application'] = p.application
		obj['customer_reference'] = p.customer_reference
		obj['project_id']='P' + p.project_id

	### 20132
	info = get_20132_info(client,project_name_swe)
	try:
                scilife_names,preps = strip_scilife_name(info.keys())
                for key in scilife_names:
                        scilife_name = scilife_names[key]
			prep = preps[key]
                        cust_name = info[key]
			incoming_QC_status = 'F' if 'F' in prep else 'P'
			try:
				obj['samples'][scilife_name] = {'customer_name': cust_name, 
								'scilife_name': scilife_name,
								'incoming_QC_status': incoming_QC_status}
			except:
				pass
        except:
                pass

	### 20158
	info = get_20158_info(client, project_name_swe)
	try:
		scilife_names, preps = strip_scilife_name(info.keys())
		duplicates = find_duplicates(scilife_names.values())
		for key in scilife_names:
			striped_scilife_name = scilife_names[key]
			status = 'inconsistent' if striped_scilife_name in duplicates else info[key][0]
			m_reads = 'inconsistent' if striped_scilife_name in duplicates else info[key][1]
			prep = preps[key]
			incoming_QC_status = 'F' if 'F' in prep else 'P'
                	if obj['samples'].has_key(striped_scilife_name):
                        	obj['samples'][striped_scilife_name]['status'] = status
                                obj['samples'][striped_scilife_name]['m_reads_sequenced'] = m_reads
			else:
				obj['samples'][striped_scilife_name] = {'scilife_name': striped_scilife_name,
									'status': status,
									'm_reads_sequenced': m_reads,
                                                               		'incoming_QC_status': incoming_QC_status}
        except:
                pass
        ### Get _id for sample_run_metrics 
        info = find_samp_from_view(samp_db, project_name)
        if len(info.keys()) > 0:
                logger.debug('sample_run_metrics found on couchdb for project %s' % project_name)
        else:
                logger.warning('No sample_run_metrics found for project %s' % project_name)
        for key in info:
                sci_name_raw = info[key][1]
                scilife_name, preps = strip_scilife_name([sci_name_raw])
                scilife_name = scilife_name[sci_name_raw]
                prep = 'A' if preps[sci_name_raw].replace('F','') == '' else preps[sci_name_raw].replace('F','')
                if obj['samples'].has_key(scilife_name):
                        if obj['samples'][scilife_name].has_key("library_prep"):
                                if obj['samples'][scilife_name]["library_prep"].has_key(prep):
                                        obj['samples'][scilife_name]["library_prep"][prep]["sample_run_metrics"][info[key][0]]=key
                                else:
                                        obj['samples'][scilife_name]["library_prep"][prep]={"sample_run_metrics":{info[key][0]:key}}
                        else:
                                obj['samples'][scilife_name]["library_prep"]={prep:{"sample_run_metrics":{info[key][0]:key}}}


	### 20135
	if WS_projects.has_key(project_name):
		logger.debug('project run on Work Set')
		info = WS_projects[project_name]	
	else:
		info={}
	info = get_20135_info(client,project_name_swe, info)
	scilife_names, preps = strip_scilife_name(info.keys())
	for key in scilife_names:
	       	striped_scilife_name = scilife_names[key]
		for prep in info[key]:
			if (obj['application']=='Exome capture')|(obj['application']=='Custom capture'):
				try: Av_sice = int(float(info[key][prep][2]))
				except: Av_sice = info[key][prep][2]
				prep_status = info[key][prep][3]
			else:
				try: Av_sice = int(float(info[key][prep][0]))
				except: Av_sice = info[key][prep][0]
				prep_status = info[key][prep][1]
                        if obj['samples'].has_key(striped_scilife_name):
                                if obj['samples'][striped_scilife_name].has_key("library_prep"):
                                        if obj['samples'][striped_scilife_name]["library_prep"].has_key(prep):
                                                obj['samples'][striped_scilife_name]["library_prep"][prep]["average_size_bp"]=Av_sice
                                                obj['samples'][striped_scilife_name]["library_prep"][prep]["prep_status"]=prep_status
                                        else:
                                                obj['samples'][striped_scilife_name]["library_prep"][prep]={"average_size_bp":Av_sice,"prep_status":prep_status}
                                else:
                                        obj['samples'][striped_scilife_name]["library_prep"]={prep:{"average_size_bp":Av_sice,"prep_status":prep_status}}
      	return obj



#		LOGGING
def my_logging(log_file):
        logger = logging.getLogger("logger")
        logger.setLevel(logging.DEBUG)

        # file handler
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        
        # console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # formatter
        formatter = logging.Formatter(
        "%(asctime)s (%(levelname)s) : %(message)s")
        fh.setFormatter(formatter)

        # add handlers to logger
        logger.addHandler(ch)
        logger.addHandler(fh)

        return logger


#		GOOGLE DOCS
def get_google_document(ssheet_title, wsheet_title, client):
	ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)
	wsheet = bcbio.google.spreadsheet.get_worksheet(client, ssheet, wsheet_title)
	content = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet)
	ss_key = bcbio.google.spreadsheet.get_key(ssheet)
	ws_key = bcbio.google.spreadsheet.get_key(wsheet)
	return content, ws_key, ss_key

def make_client(CREDENTIALS_FILE):
        credentials = bcbio.google.get_credentials({'gdocs_upload': {'gdocs_credentials': CREDENTIALS_FILE}})
        client = bcbio.google.spreadsheet.get_client(credentials)
        return client

def get_column(ssheet_content, header, col_cond=0):
	colindex=''
        for j, row in enumerate(ssheet_content):
		if colindex == '':
                	for i, col in enumerate(row):
                	        if col_cond <= i and colindex == '':
                	                if str(col).strip().replace('\n','').replace(' ','') == header.replace(' ',''):
                	                        colindex = i
		else:
                        rowindex = j-1
                        return rowindex, colindex


#		COUCHDB
def save_couchdb_obj(db, obj):
    dbobj = db.get(obj['_id'])
    time_log = datetime.utcnow().isoformat() + "Z"
    if dbobj is None:
        obj["creation_time"] = time_log 
        obj["modification_time"] = time_log 
        db.save(obj)
	return 'Created'
    else:
        obj["_rev"] = dbobj.get("_rev")
	del dbobj["modification_time"]
	obj["creation_time"] = dbobj["creation_time"]
        if not comp_obj(obj, dbobj):    
            obj["modification_time"] = time_log 
            db.save(obj)
	    return 'Uppdated'
    return None 

def comp_obj(obj, dbobj):
	for key in dbobj:
		if (obj.has_key(key)):
			if (obj[key] != dbobj[key]):
	                     return False
	     	else:
			return False
	return True

def find_proj_from_view(proj_db, project_name):
	view = proj_db.view('project/project_name')
	for proj in view:
		if proj.key == project_name:
			return proj.value
	return None

def find_samp_from_view(samp_db, proj_name):
        view = samp_db.view('names/id_to_proj')
	samps = {}
        for doc in view:
                if (doc.value[0] == proj_name)|(doc.value[0] == proj_name.lower()):
			samps[doc.key] = doc.value[1:3]
        return samps


def get_20132_info(client,project_name_swe):
        info = {}
        versions = {"01": ["Data", 'Sample name Scilife (Index included)'],
                    "02": ["Sheet1", 'Sample name Scilife'],
                    "04": ["Reception control", 'Complete sample name'],
                    "05": ["Reception control", 'SciLifeLab ID']}

        feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + '_20132', False)
        if not len(feed.entry) == 0:
                ssheet = feed.entry[0].title.text
                version = ssheet.split('_20132_')[1].split(' ')[0].split('_')[0]
                wsheet = versions[version][0]
                header = versions[version][1]
                content, ws_key, ss_key = get_google_document(ssheet, wsheet, client)
                dummy, customer_names_colindex = get_column(content, 'Sample name from customer')
                row_ind, scilife_names_colindex = get_column(content, header)
                for j,row in enumerate(content):
                        if (j > row_ind):
                                try:
                                        cust_name = str(row[customer_names_colindex]).strip()
                                        sci_name = str(row[scilife_names_colindex]).strip()
                                        if cust_name != '':
                                                info[sci_name] = cust_name
                                except:
                                        pass
                logger.debug('Google document found:    '+ssheet)
        else:
                logger.debug('Google document 20132 not found')
        return info

def get_20158_info(client, project_name_swe):
        versions = {"01": ['Sample name Scilife', "Total reads per sample", "Sheet1","Passed=P/ not passed=NP*"],
                    "02": ["Sample name (SciLifeLab)", "Total number of reads (Millions)","Sheet1",
                          "Based on total number of reads after mapping and duplicate removal"],
                    "03": ["Sample name (SciLifeLab)", "Total number of reads (Millions)","Sheet1",
                          "Based on total number of reads after mapping and duplicate removal "],
                    "05": ["Sample name (from Project read counts)", "Total number","Sheet1",
                          "Based on total number of reads","Based on total number of reads after mapping and duplicate removal"],
                    "06": ["Sample name (from Project read counts)", "Total number","Sheet1",
                          "Based on total number of reads","Based on total number of reads after mapping and duplicate removal"]}			
        info = {}
        feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + '_20158', False)
        if len(feed.entry) != 0:
                ssheet = feed.entry[0].title.text
                version = ssheet.split(str('_20158_'))[1].split(' ')[0].split('_')[0]
                content, ws_key, ss_key = get_google_document(ssheet,  versions[version][2], client)
                dummy, P_NP_colindex = get_column(content, versions[version][3])
                dummy, No_reads_sequenced_colindex = get_column(content, versions[version][1])
		row_ind, scilife_names_colindex = get_column(content, versions[version][0])
                if (version=="05")| (version=="06"):
			dummy, P_NP_duprem_colindex = get_column(content, versions[version][4]) ## [version][4] for dup rem
		else:
			P_NP_duprem_colindex=''
                for j, row in enumerate(content):
                        if (j > row_ind):
                                try:
                                        sci_name = str(row[scilife_names_colindex]).strip()
                                        no_reads = str(row[No_reads_sequenced_colindex]).strip()
                                        if (P_NP_duprem_colindex!='') and (str(row[P_NP_duprem_colindex]).strip()!=''):
                                                status = str(row[P_NP_duprem_colindex]).strip()
                                        else:
                                                status = str(row[P_NP_colindex]).strip()
                                        info[sci_name] = [status,no_reads]
                                except:
                                        pass
                logger.debug('Google document found:    '+ssheet)
        else:
                logger.debug('Google document 20158 not found')
        return info

def get_WS_info(client):
        feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, '20135', False)
	WS_projects = {}
        for ssheet in feed.entry:
                ssheet_name=ssheet.title.text
                if ssheet_name[0:3].lower()== 'ws1':
                        try:
                                content, ws_key, ss_key = get_google_document(ssheet_name, "Library QC",client)
				dummy, Transfer_to_plan_col  = get_column(content,'Transfer to planning of instrument runs')
                                dummy, Project_name_col = get_column(content,'Project name',Transfer_to_plan_col)
                                dummy, Finished_library_col  = get_column(content, 'Finished library ')
                                dummy, Hyb_Finished_library_col  = get_column(content, 'Hybridized finished library (sequence capture only)')
                                dummy, Av_sice_bp_colindex = get_column(content, 'Average size (bp) ',Finished_library_col)
                                dummy, Hyb_Av_sice_bp_colindex = get_column(content, 'Average size (bp) ',Hyb_Finished_library_col)
                                dummy, prep_status_col  = get_column(content, 'Status (P/NP)',Finished_library_col)
                                dummy, Hyb_prep_status_col  = get_column(content, 'Status (P/NP)',Hyb_Finished_library_col)
                                dummy, scilife_names_colindex = get_column(content, 'SciLifeLab ID')
                                row_ind, prep_colindex = get_column(content, 'Prep version (A, B etc)')
                                for j, row in enumerate(content):
                                        proj = str(row[Project_name_col]).strip()
					if (j > row_ind):
                                        	sci_name = str(row[scilife_names_colindex]).strip()
                                                prep = str(row[prep_colindex]).strip()
                                                if prep=='': prep='A'
                                                Av_sice_Hyb = str(row[Hyb_Av_sice_bp_colindex]).strip()
                                                prep_status_Hyb = str(row[Hyb_prep_status_col]).strip()
                                                Av_sice = str(row[Av_sice_bp_colindex]).strip()
                                                prep_status = str(row[prep_status_col]).strip()
						prep_info_list = [Av_sice,prep_status,Av_sice_Hyb,prep_status_Hyb]
                                                if ''.join(prep_info_list) !='' and not WS_projects.has_key(proj):
                                                        WS_projects[proj]={sci_name:{prep:prep_info_list}}
                                                elif WS_projects[proj].has_key(sci_name) and ''.join(prep_info_list) !='':
                                                        WS_projects[proj][sci_name][prep] = prep_info_list
                                                elif ''.join(prep_info_list) !='':
                                                        WS_projects[proj][sci_name] = {prep:prep_info_list}
                        except:
                                pass
        return WS_projects

def get_20135_info(client,project_name_swe, info={}):
        feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + '_20135', False)
        if len(feed.entry) != 0:
                ssheet = feed.entry[0].title.text
                content, ws_key, ss_key = get_google_document(ssheet, "Library QC", client)
		dummy, Hyb_Finished_library_col  = get_column(content, 'Hybridized finished library (sequence capture only)')
                dummy, Finished_library_col  = get_column(content, 'Finished library ')
                dummy, Av_sice_bp_colindex = get_column(content,'Average size (bp)',Finished_library_col)
                dummy, Hyb_Av_sice_bp_colindex = get_column(content, 'Average size (bp) ',Hyb_Finished_library_col)
                dummy, prep_status_col  = get_column(content, 'Status (P/NP)',Finished_library_col)
                dummy, Hyb_prep_status_col  = get_column(content, 'Status (P/NP)',Hyb_Finished_library_col)
                row_ind, scilife_names_colindex = get_column(content,'SciLifeLab ID')
                row_ind, prep_colindex = get_column(content,'Prep version (A, B etc)')
                for j, row in enumerate(content):
                        if (j > row_ind):
                                try:
                                        sci_name = str(row[scilife_names_colindex]).strip()
                                        Av_sice_Hyb = str(row[Hyb_Av_sice_bp_colindex]).strip()
                                        prep_status_Hyb = str(row[Hyb_prep_status_col]).strip()
                                        Av_sice = str(row[Av_sice_bp_colindex]).strip()
                                        prep_status = str(row[prep_status_col]).strip()
                                        prep = str(row[prep_colindex]).strip()
                                        if prep=='': prep='A'
					prep_info_list = [Av_sice,prep_status,Av_sice_Hyb,prep_status_Hyb]
                                        if info.has_key(sci_name) and ''.join(prep_info_list) !='':
                                                info[sci_name][prep] = prep_info_list
                                        elif ''.join(prep_info_list) !='':
                                                info[sci_name] = {prep:prep_info_list}
                                except:
                                        pass
               	logger.debug('Google document found:	'+ssheet)
        else:
               	logger.debug('Google document 20135 not found')
	return info


#		NAME HANDELING
def _replace_ascii(str):
    # Substitute swedish characters for sensible counterparts
    str = str.replace(u'\xc5','A')
    str = str.replace(u'\xe5','a')
    str = str.replace(u'\xc4','A')
    str = str.replace(u'\xe4','a')
    str = str.replace(u'\xd6','O')
    str = str.replace(u'\xf6','o')
    return str.encode('ascii','replace')

def find_duplicates(list):
	dup = []
        shown = []
        for name in list:
        	if name in shown and name not in dup:
                	dup.append(name)
                shown.append(name)
	return dup

def strip_scilife_name(names):
	N = {}
	P = {}
        preps = 'F_BCDE'
        for name_init in names:
		prep = ''
		indexes = ['_index','_rpi','_agilent','_mondrian','_haloht','_halo','_sureselect','_dual','_hht','_ss','_i','_r','_a','_m','_h']
		name = name_init.replace('-', '_').replace(' ', '')
		for i in indexes:
			name=name.split(i)[0]
		if name != '':
			while name[-1] in preps:
				prep = name[-1] + prep
				name = name[0: -1]
			if name != '':
                        	N[name_init] = name
				P[name_init] = prep.replace('_', '')
	return N, P


def  main(client, CONFIG, URL, proj_ID, all_projects, GPL ):
	couch = couchdb.Server("http://" + URL)
       	samp_db = couch['samples']
        proj_db = couch['projects']
	info = None
	WS_projects = get_WS_info(client)
        if all_projects:
		content, ws_key, ss_key = get_google_document("Genomics Project list", GPL, client) 	
		row_ind, col_ind = get_column(content, 'Project name')
		for j, row in enumerate(content):
      			try:
	        		proj_ID = str(row[col_ind]).strip().split(' ')[0]
				if (proj_ID != '') & (j > row_ind + 2):
                       			obj = get_proj_inf(WS_projects,proj_ID, samp_db, proj_db, client, CONFIG)
        				if obj['samples'].keys() != []:
                				info = save_couchdb_obj(proj_db, obj)
						if info:
							logger.info('CouchDB: %s %s %s' % (obj['_id'], obj['project_name'], info))
			except:
				pass		
	elif proj_ID is not None:
	        obj = get_proj_inf(WS_projects,proj_ID, samp_db, proj_db, client, CONFIG)
        	if obj['samples'].keys() != []:
                	info = save_couchdb_obj(proj_db, obj)
	else:
		logger.debug('Argument error')
	if info:
		logger.info('CouchDB: %s %s %s' % (obj['_id'], obj['project_name'], info))

if __name__ == '__main__':
    	usage = """Usage:	python project_summary_upload.py [options]

Options (Only one option is acceptab):
	-a,			upploads all projects in genomics project list into couchDB
     	-p <project_ID>,	upploads the project <project_ID> into couchDB                                         
	-g <GPL_sheet>,		specify sheet in GPL if not Ongoing (default)
"""
    	parser = OptionParser(usage=usage)
    	parser.add_option("-p", "--project", dest="project_ID", default=None)
    	parser.add_option("-a", "--all_projects", dest="all_projects", action="store_true", default=False)
	parser.add_option("-g", "--GPL", dest="GPL", default="Ongoing")
    	(options, args) = parser.parse_args()

	CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')
	CONFIG_FILE = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')
	CONFIG = cl.load_config(CONFIG_FILE)
	log_path = CONFIG['analysis']['log']
	URL = CONFIG['couch_db']['maggie_login']+':'+CONFIG['couch_db']['maggie_pass']+'@'+CONFIG['couch_db']['maggie_url'] + ':' + str(CONFIG['couch_db']['maggie_port'])
	logger = my_logging(log_path+'/proj_coucdb.log')

	if (options.project_ID is None) and (options.all_projects is False):
                sys.exit()
	kwargs = {'proj_ID': options.project_ID,
                  'all_projects': options.all_projects,
		   'GPL':options.GPL}
    	main(make_client(CREDENTIALS_FILE), CONFIG, URL , **kwargs)



