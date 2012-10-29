#!/usr/bin/env python
"""ProjectSummaryUploadV2
Created by Maya Brandi on 2012-09-00.

Creates and updates project_summary documents in the projects database 
on couchdb. Finds projects by looping through the 20132 documents on 
google docs.

Sources of the fields in project_summary:
google docs:
	Document 20132: 
		project_id, 
		scilife_name (samples), 
		customer_name (samples) 
	Genomics project list: 
		min_m_reads_per_sample_ordered, 
		application, 
		customer_reference, 
		uppnex_id, 
		no_of_samples
	Document 20158:
		status (samples), 
		m_reads_sequenced (samples)
                scilife_name (samples), 
                customer_name (samples)
couchdb:
	samples:
		sample_run_metrics

"""
import sys
import os
import time
from  datetime  import  datetime
from uuid import uuid4
import hashlib
from optparse import OptionParser
import logging

import bcbio
import bcbio.google
import bcbio.scilifelab.google.project_metadata as pm
import bcbio.pipeline.config_loader as cl
from bcbio.google import _to_unicode, spreadsheet
import couchdb
#from scilifelab.utils.string import replace_ascii

def get_proj_inf(project_name_swe, samp_db, proj_db, credentials_file, config_file):

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
             'project_id': project_name, 
             '_id': key}


	### Get minimal #M reads and uppnexid from Genomics Project list
	logger.debug('Getting minimal #M reads and uppnexid from Genomics Project list for project %s' % project_name_swe)

	config = cl.load_config(config_file)
	p = pm.ProjectMetaData(project_name, config)
	if p.project_name is None:
		p = pm.ProjectMetaData(project_name_swe, config)
	if p.project_name is None:
		logger.warning('Google Document Genomics Project list: %s not found' % project_name) 
	else:
		if p.min_reads_per_sample.strip() != '':
                	obj['min_m_reads_per_sample_ordered'] = float(p.min_reads_per_sample)
		if p.no_samples.strip() != '':
			obj['no_of_samples'] = int(p.no_samples)
                obj['uppnex_id'] = p.uppnex_id
		obj['application'] = p.application
		obj['customer_reference'] = p.customer_reference


	### 20132
	logger.debug('Trying to find Scilife Sample names from table 20132')

       	versions = {"01": ["Data", 'Sample name Scilife (Index included)'],
		    "02": ["Sheet1", 'Sample name Scilife'],
		    "04": ["Reception control", 'Complete sample name'],
		    "05": ["Reception control", 'SciLifeLab ID']}

	# Load google document
	client = make_client(credentials_file)
	feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + '_20132', False)
	if len(feed.entry) == 0:
    		ssheet = None
		logger.warning("Could not find spreadsheet 20132 for %s" % project_name_swe)
	else:
    		ssheet = feed.entry[0].title.text
  		version	= ssheet.split('_20132_')[1].split(' ')[0].split('_')[0]
		wsheet = versions[version][0]
		header = versions[version][1]
		content, ws_key, ss_key = get_google_document(ssheet, wsheet, credentials_file)
		logger.debug("Document found")
		logger.debug(ssheet)	

	# Get Scilife Sample names
	try:    
	   	dummy, customer_names_colindex = get_column(content, 'Sample name from customer')
		row_ind, scilife_names_colindex = get_column(content, header)
		info = {}
                for j,row in enumerate(content):
			if (j > row_ind):
				try:
					cust_name = str(row[customer_names_colindex]).strip()
					sci_name = str(row[scilife_names_colindex]).strip()	
					if cust_name != '':
						info[sci_name] = cust_name
				except:
					pass

		logger.debug('Names found')
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
		logger.debug('Names not found')
                pass

	### 20158
	logger.debug('Getting Sample Status from table 20158')

        versions = {"01": ['Sample name Scilife', "Total reads per sample", 
			"Passed=P/ not passed=NP*", 'Sample name from customer'],
                    "02": ["Sample name (SciLifeLab)", "Total number of reads (Millions)", 
			"Based on total number of reads", 'Sample name (customer)'],
                    "03": ["Sample name (SciLifeLab)", "Total number of reads (Millions)", 
			"Based on total number of reads", 'Sample name (customer)']}

        # Load google document
	mistakes = ["_", " _", " ", ""]
	found = FALSE
	for m in mistakes:
		feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + m + '20158', False)
        	if len(feed.entry) != 0:
			try:
				ssheet = feed.entry[0].title.text
				version = ssheet.split(str(m + '20158_'))[1].split(' ')[0].split('_')[0]	
				content, ws_key, ss_key = get_google_document(ssheet, "Sheet1", credentials_file)
				found = TRUE
				break
                	except:
				pass
	if found:
		logger.debug('Google document found')
		logger.debug(ssheet)
	else:
		logger.warning("Could not find spreadsheet 20158 for %s" % project_name_swe)

	# Get status etc from loaded document
	try:
		dummy, P_NP_colindex = get_column(content, versions[version][2])
		dummy, No_reads_sequenced_colindex = get_column(content, versions[version][1])
		dummy, customer_names_colindex = get_column(content, versions[version][3])
        	row_ind, scilife_names_colindex = get_column(content, versions[version][0])
		info = {}
                for j, row in enumerate(content):
			if (j > row_ind):
				try:
                                        sci_name = str(row[scilife_names_colindex]).strip()
					cust_name = str(row[customer_names_colindex]).strip()
                                        no_reads = str(row[No_reads_sequenced_colindex]).strip()
                                        status = str(row[P_NP_colindex]).strip()
                                        info[sci_name] = [status,no_reads,cust_name]
				except:
					pass
		scilife_names, preps = strip_scilife_name(info.keys())
		duplicates = find_duplicates(scilife_names.values())
		for key in scilife_names:
			striped_scilife_name = scilife_names[key]
			status = 'inconsistent' if striped_scilife_name in duplicates else info[key][0]
			m_reads = 'inconsistent' if striped_scilife_name in duplicates else info[key][1]
			cust_name = info[key][2]
			prep = preps[key]
			incoming_QC_status = 'F' if 'F' in prep else 'P'
                	if obj['samples'].has_key(striped_scilife_name):
                        	obj['samples'][striped_scilife_name]['status'] = status
                                obj['samples'][striped_scilife_name]['m_reads_sequenced'] = m_reads
			else:
				obj['samples'][striped_scilife_name] = {'customer_name': cust_name, 
									'scilife_name': striped_scilife_name,
									'status': status,
									'm_reads_sequenced': m_reads,
                                                               		'incoming_QC_status': incoming_QC_status}
			
        except:
                pass


        ### Get _id for sample_run_metrics 
        logger.debug('Getting _id for sample_run_metrics')
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
        logger.debug('Getting average read length from table 20135')

	versions = {"04": ['SciLifeLab ID','Prep version (A, B etc)','Average size (bp)'],
	            "05": ['SciLifeLab ID','Prep version (A, B etc)','Average size (bp)'],
	            "06": ['SciLifeLab ID','Prep version (A, B etc)','Average size (bp)']}
	mistakes = ["_","_ ", " _", " ", ""]
	found = FALSE

	for m in mistakes:
	        feed = bcbio.google.spreadsheet.get_spreadsheets_feed(client, project_name_swe + m + '20135', False)
	        if len(feed.entry) != 0:
	                ssheet = feed.entry[0].title.text
			version = ssheet.split('20135')[1].replace('_',' ').lstrip(' ').split(' ')[0]
	                content, ws_key, ss_key = get_google_document(ssheet, "Library QC", credentials_file)
	                found = TRUE
                
	if found:
        	logger.debug('Google document found')
		logger.debug(ssheet)
	else:
	        logger.debug('Google document not found')

	# Get average read length from loaded document
	try:
		dummy, Finished_library_col  = get_column(content, 'Finished library ')
		dummy, Av_sice_bp_colindex = get_column(content, versions[version][2],Finished_library_col)
		row_ind, scilife_names_colindex = get_column(content, versions[version][0])
		row_ind, prep_colindex = get_column(content, versions[version][1])
		info = {}
		for j, row in enumerate(content):
		        if (j > row_ind):
		                try:
		                        sci_name = str(row[scilife_names_colindex]).strip()
		                        Av_sice = str(row[Av_sice_bp_colindex]).strip()
		                        prep = str(row[prep_colindex]).strip()
		                        info[sci_name] = [Av_sice,prep]
		                except:
		                        pass
		scilife_names, preps = strip_scilife_name(info.keys())
	        for key in scilife_names:
	        	striped_scilife_name = scilife_names[key]
	                Av_sice = info[key][0]
			if info[key][1].strip() != '':
				prep=info[key][1]#    KOntrollera!!!!!
			elif preps[key].strip()!= '':
				prep=preps[key]
				prep = 'A' if preps[key].replace('F','') == '' else preps[key].replace('F','')
	                try:
	               		obj['samples'][striped_scilife_name]["library_prep"][prep]["average_size_bp"]=Av_sice
        	        except:
	                	pass
        except:
                pass
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
def get_google_document(ssheet_title, wsheet_title, credentials_file):
	client = make_client(credentials_file)
	ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)
	wsheet = bcbio.google.spreadsheet.get_worksheet(client, ssheet, wsheet_title)
	content = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet)
	ss_key = bcbio.google.spreadsheet.get_key(ssheet)
	ws_key = bcbio.google.spreadsheet.get_key(wsheet)
	return content, ws_key, ss_key

def make_client(credentials_file):
        credentials = bcbio.google.get_credentials({'gdocs_upload': {'gdocs_credentials': credentials_file}})
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

def find_proj_from_view(proj_db, proj_id):
	view = proj_db.view('project/project_id')
	for proj in view:
		if proj.key == proj_id:
			return proj.value
	return None

def find_samp_from_view(samp_db, proj_id):
        view = samp_db.view('names/id_to_proj')
	samps = {}
        for doc in view:
                if doc.value[0] == proj_id:
			samps[doc.key] = doc.value[1:3]
        return samps

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
		name = name_init.replace('-', '_').replace(' ', '').split("_index")[0].strip()
		if name != '':
			while name[-1] in preps:
				prep = name[-1] + prep
				name = name[0: -1]
			if name != '':
                        	N[name_init] = name
				P[name_init] = prep.replace('_', '')
	return N, P


def  main(credentials_file, config_file, URL, proj_ID, all_projects):
	client = make_client(credentials_file)
	couch = couchdb.Server("http://" + URL)
       	samp_db = couch['samples']
        proj_db = couch['projects']
	info = None
	
        if all_projects:
		content, ws_key, ss_key = get_google_document("Genomics Project list", "Ongoing", credentials_file)	
		row_ind, col_ind = get_column(content, 'Project name')
		for j, row in enumerate(content):
      			try:
	        		proj_ID = str(row[col_ind]).strip().split(' ')[0]
				if (proj_ID != '') & (j > row_ind + 2):
                       			obj = get_proj_inf(proj_ID, samp_db, proj_db, credentials_file, config_file)
        				if obj['samples'].keys() != []:
                				info = save_couchdb_obj(proj_db, obj)
						if info:
							logger.debug('couchdb %s' % info)
			except:
				pass		
	elif proj_ID is not None:
	        obj = get_proj_inf(proj_ID, samp_db, proj_db, credentials_file, config_file)
        	if obj['samples'].keys() != []:
                	info = save_couchdb_obj(proj_db, obj)
	else:
		logger.debug('Argument error')
	if info:
		logger.debug('couchdb %s' % info)
		logger.info('CouchDB: %s %s %s' % (obj['_id'], obj['project_id'], info))

if __name__ == '__main__':
    	usage = """Usage:	python project_summary_upload.py <url> [options]

Arguments:
	<url>			Database url (excluding http://)
Options (Only one option is acceptab):
	-a,			upploads all projects in genomics project list into couchDB
     	-p <project_ID>,	upploads the project <project_ID> into couchDB                                         

"""
	logger = my_logging('proj_coucdb.log')
    	parser = OptionParser(usage=usage)
    	parser.add_option("-p", "--project", dest="project_ID", default=None)
    	parser.add_option("-a", "--all_projects", dest="all_projects", action="store_true", default=False)
    	(options, args) = parser.parse_args()

	CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')
	CONFIG_FILE = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')


	if (options.project_ID is None) and (options.all_projects is False):
                print usage
                sys.exit()
	kwargs = {'proj_ID': options.project_ID,
                  'all_projects': options.all_projects}

    	main(CREDENTIALS_FILE, CONFIG_FILE, *args, **kwargs)



