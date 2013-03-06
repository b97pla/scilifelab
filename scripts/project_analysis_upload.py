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

def get_analysis_inf(project_name_swe, analysis_db, proj_db,config):
	project_name = _replace_ascii(_to_unicode(project_name_swe))

	key = find_proj_from_view(analysis_db, project_name)
	if not key: key = uuid4().hex
        logger.info('Handling proj %s %s' % (project_name, key))

        obj={'application': '',
             'no_of_samples': '',
             'entity_type': 'project_analysis',
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
                if p.no_samples.strip() != '':
                        obj['no_of_samples'] = int(p.no_samples)
                obj['application'] = p.application
                obj['project_id']='P' + p.project_id

	### Files from analysis
	stat = open('stat.json','r')
	RSeQC = open('RSeQC_rd.json','r')
	stat_dict = eval(stat.read())
	RSeQC_dict = eval(RSeQC.read())
	scilife_names,preps = strip_scilife_name(RSeQC_dict.keys())
	for samp in scilife_names:
		samp_name_stripped = scilife_names[samp]
		prep = preps[samp]
		picardDup = pars_picardDup_metrics('tophat_out_' + samp + '/' + samp + '_picardDup_metrics')
		tophat_version = pars_tophat_reports_log('tophat_out_' + samp + '/logs/reports.log' )
		obj['samples'][samp_name_stripped] = 	{'mapping_statistics': stat_dict[samp],
							'read_distribution': RSeQC_dict[samp],
							'picard_dup': picardDup,
							'prep':prep,
							'tophat_version':tophat_version}
	
        rRNA = open('rRNA.quantification','r')
        for line in rRNA: 
		samp 		    = line.split()[0]
		scilife_names,preps = strip_scilife_name([samp])
		samp_name_stripped  = scilife_names[samp]
		obj['samples'][samp_name_stripped]['percent_rRNA'] = float(line.split()[1].strip('%'))
      	return obj

def pars_tophat_reports_log(file):
	f=open(file,'r')
	for line in f:
		if line.split()[0]=='tophat_reports':
			return line.split()[1]
def pars_picardDup_metrics(file):
	f=open(file,'r')
	lines=f.readlines()
	for i,line in enumerate(lines):
		if len(line.split('\t'))>0 and line.split('\t')[0]=='## METRICS CLASS':
			keys=lines[i+1].strip().split('\t')
		        vals=lines[i+2].strip().split('\t')
	                return dict(zip(keys,vals))

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

def strip_scilife_name(names):
	N = {}
	P = {}
        preps = 'F_BCDE'
        for name_init in names:
		prep = ''
		name = name_init.replace('-', '_').replace(' ', '').split("_index")[0].split("_ss")[0].split("_dual")[0].strip()
		if name != '':
			while name[-1] in preps:
				prep = name[-1] + prep
				name = name[0: -1]
			if name != '':
                        	N[name_init] = name
				P[name_init] = prep.replace('_', '')
	return N, P


def  main(client, config, URL, proj_ID):
	couch = couchdb.Server("http://" + URL)
       	analysis_db = couch['analysis']
        proj_db = couch['projects']
	
	obj = get_analysis_inf(proj_ID, analysis_db, proj_db,config)
        if obj['samples'].keys() != []:
                info = save_couchdb_obj(analysis_db, obj)
		try:
			logger.info('CouchDB: %s %s %s' % (obj['_id'], obj['project_name'], info))
		except:
			pass

if __name__ == '__main__':
    	usage = """Usage:	python project_analysis_upload.py <project_ID> """

	CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')
	client = make_client(CREDENTIALS_FILE)
	CONFIG_FILE = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')
	CONFIG = cl.load_config(CONFIG_FILE)
	log_path = CONFIG['analysis']['log']
	URL = CONFIG['couch_db']['maggie_url'] + ':' + str(CONFIG['couch_db']['maggie_port'])
	logger = my_logging(log_path+'/analysis_coucdb.log')

	if len(sys.argv) != 2:
                sys.exit(usage)
	roject_ID=sys.argv[1]
    	main(client, CONFIG, URL , roject_ID)

