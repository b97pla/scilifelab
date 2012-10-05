#!/usr/bin/env python
"""
ProjectSummaryUploadV2
Created by Maya Brandi on 2012-09-00.

Creates and updates project_summary documents in the projects database on couchdb. 
Finds projects by looping through the 20132 documents on google docs.

Sources of the fields in project_summary:
google docs:
	Document 20132: 
		project_id, scilife_name (samples), customer_name (samples) 
	Genomics project list: 
		min_m_reads_per_sample_ordered, application, customer_reference, uppnex_id, no_of_samples
	Document 20158:
		status (samples), m_reads_sequenced (samples)
couchdb:
	samples:
		sample_run_metrics
"""
from uuid import uuid4
import bcbio.google
import bcbio.google.spreadsheet
from optparse import OptionParser
import sys
import os
import hashlib
import couchdb
import time
from  datetime  import  datetime
#from scilifelab.scripts.process_run_info import _replace_ascii
import bcbio.scilifelab.google.project_metadata as pm
import bcbio.pipeline.config_loader as cl
import logging
from bcbio.google import _to_unicode

def get_proj_inf(project_name_swe,samp_db,proj_db,credentials_file,config_file):
	logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',filename='proj_coucdb.log',level=logging.INFO)
	project_name	 = _replace_ascii(_to_unicode(project_name_swe))
	key		= find_proj_from_view(proj_db,project_name)
	if not key:
		key = uuid4().hex

        logging.info(str('Handling proj '+project_name+' '+ key))
        print key

        obj={   'application':'',
		'customer_reference':'',
		'min_m_reads_per_sample_ordered':'',
		'no_of_samples':'',
                'entity_type': 'project_summary',
                'uppnex_id': '',                 
		'samples': {},
                'project_id': project_name, 
                '_id': key}


	### Get minimal #M reads and uppnexid from Genomics Project list
	print '\nGetting minimal #M reads and uppnexid from Genomics Project list for project ' + project_name_swe

	config = cl.load_config(config_file)
	p = pm.ProjectMetaData(project_name,config)
	if p.project_name == None:
		p = pm.ProjectMetaData(project_name_swe,config)
	if p.project_name == None:
		print project_name+' not found in genomics project list'
		logging.warning(str('Google Document Genomics Project list: '+project_name+' not found')) 
	else:
		if p.min_reads_per_sample.strip() !='':
                	obj['min_m_reads_per_sample_ordered'] = float(p.min_reads_per_sample)
		if p.no_samples.strip() !='':
			obj['no_of_samples'] = int(p.no_samples)
                obj['uppnex_id']                      = p.uppnex_id
		obj['application']		      = p.application
		obj['customer_reference']             = p.customer_reference


	### Get costumer and Scilife Sample name from _20132_0X_Table for Sample Summary and Reception Control
	print '\nTrying to find Scilife Sample names from '+project_name_swe+'_20132_0X_Table for Sample Summary and Reception Control'

       	versions = { 	"01":["Data",'Sample name Scilife (Index included)'],
			"02":["Sheet1",'Sample name Scilife'],
			"04":["Reception control",'Complete sample name'],
			"05":["Reception control",'SciLifeLab ID']}

	# Load google document
	client	= make_client(credentials_file)
	feed 	= bcbio.google.spreadsheet.get_spreadsheets_feed(client,project_name_swe+'_20132', False) #FIXA: Hantera mistakes
	if len(feed.entry) == 0:
    		ssheet=None
		logging.warning("Google Document %s: Could not find spreadsheet" % str(project_name_swe+'_20132_XXX'))
		print "Could not find spreadsheet" 
	else:
    		ssheet	= feed.entry[0].title.text
  		version	= ssheet.split('_20132_')[1].split(' ')[0].split('_')[0]
		wsheet 	= versions[version][0]
		header 	= versions[version][1]
		content, ws_key, ss_key = get_google_document(ssheet, wsheet, credentials_file)
	
	# Get Scilife Sample names
	try:    
	   	dummy, customer_names_colindex 	= get_column(content,'Sample name from customer')
		row_ind, scilife_names_colindex = get_column(content, header)
		info={}
                for j,row in enumerate(content):
			if (j > row_ind):
				try:
					cust_name = str(row[customer_names_colindex]).strip()
					sci_name  = str(row[scilife_names_colindex]).strip().replace('-','_')
					if cust_name != '':
						info[sci_name] = cust_name
				except:
					pass
		print 'Names found'
		for scilife_name in info:
			try:
				obj['samples'][scilife_name] = {'customer_name': info[scilife_name], 'scilife_name':scilife_name}
			except:
				pass
        except:
		print 'Names not found'
                pass

	### Get Sample Status from _20158_01_Table for QA HiSeq2000 sequencing results for samples
	print '\nGetting Sample Status from '+project_name_swe+'_20158_0X_Table for QA HiSeq2000 sequencing results for samples'

        versions = {    "01":['Sample name Scilife',"Total reads per sample","Passed=P/ not passed=NP*",'Sample name from customer'],
                        "02":["Sample name (SciLifeLab)","Total number of reads (Millions)","Based on total number of reads",'Sample name (customer)'],
                        "03":["Sample name (SciLifeLab)","Total number of reads (Millions)","Based on total number of reads",'Sample name (customer)']}

        # Load google document
	mistakes = ["_"," _"," ",""]
	found='FALSE'
	for m in mistakes:
		feed    = bcbio.google.spreadsheet.get_spreadsheets_feed(client,project_name_swe + m + '20158', False)
        	if len(feed.entry) != 0:
			try:
				ssheet  = feed.entry[0].title.text
				version = ssheet.split(str(m+'20158_'))[1].split(' ')[0].split('_')[0]	
				content, ws_key, ss_key = get_google_document(ssheet,"Sheet1",credentials_file)
				found='TRUE'
				break
                	except:
				pass
	if found=='TRUE':
		print 'Google document found!'
	else:
		print 'Google document NOT found!'
		logging.warning("Google Document %s: Could not find spreadsheet" % str(project_name_swe+'_20158_XXX'))

	# Get status etc from loaded document
	try:
		dummy, P_NP_colindex 			= get_column(content,versions[version][2])
		dummy, No_reads_sequenced_colindex 	= get_column(content,versions[version][1])
		dummy, customer_names_colindex          = get_column(content,versions[version][3])
        	row_ind, scilife_names_colindex 	= get_column(content,versions[version][0])
		info={}
                for j,row in enumerate(content):
			if ( j > row_ind ):
				try:
                                        sci_name=str(row[scilife_names_colindex]).strip()
					cust_name=str(row[customer_names_colindex]).strip()
                                        no_reads=str(row[No_reads_sequenced_colindex]).strip()
                                        if sci_name[-1]=='F':
                                                status='P'
                                        else:
                                                status	=str(row[P_NP_colindex]).strip()
                                        info[sci_name]	=[status,no_reads,cust_name]
				except:
					pass
		scilife_names 	= strip_scilife_name(info.keys())
		duplicates	= find_duplicates(scilife_names.values())
		for key in scilife_names:
			striped_scilife_name 	= scilife_names[key]
			status			= info[key][0]
			m_reads			= info[key][1]
			cust_name		= info[key][2]
			if striped_scilife_name in duplicates:
                                        status 	= 'inconsistent'
					m_reads	= 'inconsistent'
			try:
                		if obj['samples'].has_key(striped_scilife_name):
                        		obj['samples'][striped_scilife_name]['status']            = status
                        	        obj['samples'][striped_scilife_name]['m_reads_sequenced'] = m_reads
				else:
					obj['samples'][striped_scilife_name]=	{'customer_name': cust_name, 
										'scilife_name':striped_scilife_name,
										'status':status,
										'm_reads_sequenced':m_reads}
			except:
				pass
        except:
		print 'Status and M reads sequenced not found in '+project_name_swe+'_20158_0X_Table for QA HiSeq2000 sequencing results for samples'
                pass


	### Get _id for sample_run_metrics and bcbb names -- use couchdb views instead.... To be fixed...
	print '\nGetting _id for sample_run_metrics'

	info	= find_samp_from_view(samp_db,project_name)

	if len(info.keys())>0:
		print 'sample_run_metrics found on couchdb for project '+ project_name
	else:
		print 'no sample_run_metrics found on couchdb for project '+ project_name
		logging.warning(str('CouchDB: No sample_run_metrics found for project '+ project_name))
        for key in info:
        	scilife_name = strip_scilife_name([info[key][1]])[info[key][1]]
                if obj['samples'].has_key(scilife_name):
        		if obj['samples'][scilife_name].has_key("sample_run_metrics"):
                		obj['samples'][scilife_name]["sample_run_metrics"][info[key][0]]=key
                        else:
                              	obj['samples'][scilife_name]["sample_run_metrics"] = {info[key][0]:key}
	return obj

#		GOOGLE DOCS
def get_google_document(ssheet_title,wsheet_title,credentials_file):
	""""""
	client  = make_client(credentials_file)
	ssheet 	= bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
	wsheet 	= bcbio.google.spreadsheet.get_worksheet(client,ssheet,wsheet_title)
	content = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet)
	ss_key 	= bcbio.google.spreadsheet.get_key(ssheet)
	ws_key 	= bcbio.google.spreadsheet.get_key(wsheet)
	return content, ws_key, ss_key

def make_client(credentials_file):
        credentials = bcbio.google.get_credentials({'gdocs_upload': {'gdocs_credentials': credentials_file}})
        client  = bcbio.google.spreadsheet.get_client(credentials)
        return client

def get_column(ssheet_content,header):
	""""""
	colindex=''
	for j,row in enumerate(ssheet_content):
                if colindex == '':
			try:
                        	for i,col in enumerate(row):
					try:
                                		if str(col).strip() == header:
                                        		colindex = i
					except:
						pass
			except:
				pass
		else:
			rowindex = j-1
			return rowindex, colindex

#		COUCHDB
def save_couchdb_obj(db, obj):
    dbobj	= db.get(obj['_id'])
    time_log 	= datetime.utcnow().isoformat() + "Z"
    if dbobj is None:
        obj["creation_time"] 	 = time_log 
        obj["modification_time"] = time_log 
        db.save(obj)
	return 'Created'
    else:
        obj["_rev"] = dbobj.get("_rev")
	del dbobj["modification_time"]
	obj["creation_time"] = dbobj["creation_time"]
        if comp_obj(obj,dbobj)==False:    
            obj["modification_time"] = time_log 
            db.save(obj)
	    return 'Uppdated'
    return None 

def comp_obj(obj,dbobj):
	for key in dbobj:
		if (obj.has_key(key)):
			if (obj[key]!=dbobj[key]):
	                     return False
	     	else:
			return False
	return True

def find_proj_from_view(proj_db,proj_id):
	view = proj_db.view('project/project_id')
	for proj in view:
		if proj.key==proj_id:
			return proj.value
	return None

def find_samp_from_view(samp_db,proj_id):
        view = samp_db.view('names/id_to_proj')
	samps={}
        for doc in view:
                if doc.value[0]==proj_id:
			samps[doc.key]=doc.value[1:3]
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
	dup=[]
        shown=[]
        for name in list:
        	if name in shown and name not in dup:
                	dup.append(name)
                shown.append(name)
	return dup

def strip_scilife_name(names):
	N={}
        preps = '_BCDE'
        for name_init in names:
		name = name_init.replace('-','_').split(' ')[0].split("_index")[0].strip()
		if name !='':
			while name[-1] in preps:
				name=name[0:-1]
			if name !='':
                        	N[name_init]=name
	return N


def  main(credentials_file,config_file, proj_ID,all_projects):
	client	= make_client(credentials_file)
	couch   = couchdb.Server("http://maggie.scilifelab.se:5984")
       	samp_db = couch['samples']
        proj_db = couch['projects']
	info	= None
	
        if all_projects == True:
		content, ws_key, ss_key = get_google_document("Genomics Project list","Ongoing",credentials_file)	
		row_ind, col_ind = get_column(content,'Project name')
		for j,row in enumerate(content):
      			try:
	        		proj_ID = str(row[col_ind]).strip().split(' ')[0]
				if (proj_ID !='') & (j>row_ind+2):
                       			obj     = get_proj_inf(proj_ID ,samp_db,proj_db ,credentials_file, config_file)
        				if obj['samples'].keys()!=[]:
                				info =  save_couchdb_obj(proj_db, obj)
			except:
				pass		
	elif proj_ID !=None:
	        obj = get_proj_inf(proj_ID ,samp_db ,proj_db ,credentials_file, config_file )
        	if obj['samples'].keys()!=[]:
                	info = save_couchdb_obj(proj_db, obj)
	else:
		print 'Argument error'
	if info:
		print 'couchdb '+info
        	logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',filename='proj_coucdb.log',level=logging.INFO)
		logging.info('CouchDB: '+obj['_id'] + ' ' + obj['project_id'] + ' ' +info)


if __name__ == '__main__':
    	usage = """
Usage:	python ProjectSummaryUploadV2.py [options]

One and only one option is acceptable:
	-a,			upploads all projects in genomics project list into couchDB
     	-p <project_ID>,	upploads the project <project_ID> into couchDB                                         
"""

    	parser = OptionParser(usage=usage)
    	parser.add_option("-p", "--project", dest="project_ID", default=None)
    	parser.add_option("-a", "--all_projects", dest="all_projects",action="store_true", default=False)
    	(options, args) = parser.parse_args()

	CREDENTIALS_FILE = os.path.join(os.environ['HOME'],'opt/config/gdocs_credentials')
	CONFIG_FILE	 = os.path.join(os.environ['HOME'],'opt/config/post_process.yaml')


	if (options.project_ID==None) and (options.all_projects==False):
                print usage
                sys.exit()

    	main(CREDENTIALS_FILE,CONFIG_FILE,options.project_ID,options.all_projects)















