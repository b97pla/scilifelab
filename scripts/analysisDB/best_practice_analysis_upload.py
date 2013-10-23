#!/usr/bin/env python

"""Script to load runinfo from the lims process: 'Illumina Sequencing (Illumina SBS) 4.0' 
into the flowcell database in statusdb.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""
import sys
import os
import codecs
from optparse import OptionParser
from pprint import pprint
from datetime import date
from scilifelab.db.statusDB_utils import *
import scilifelab.log
import analysisDB as DB

from scilifelab.google.google_docs import SpreadSheet
from scilifelab.google import get_credentials

def main(project_name, conf, cred):
    credentials = get_credentials(cred)
    client = SpreadSheet(credentials)
    config = cl.load_config(conf)
    couch = load_couch_server(conf)
    analysis_db = couch['analysis']
    print dir(analysis_db)
    #proj_db = couch['analysis']#['projects']
    BP_RNA = DB.BP_RNA(project_name)
    key = find_proj_from_view(analysis_db, project_name)
    BP_RNA.obj['_id'] = find_or_make_key(key)
    info = save_couchdb_obj(analysis_db, BP_RNA.obj)
    LOG.info('project %s %s : _id = %s' % (project_name, info, BP_RNA.obj['_id']))

if __name__ == '__main__':
    usage = "Usage:       python best_practice_analysis_upload.py [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--project", dest="project_name", default=None, 
    help = "eg: J.Doe_13_01.")
    
    parser.add_option("-g", "--credentials", dest="cred",
    default=os.path.join(os.environ['HOME'],'opt/config/gdocs_credentials'),
    help = "Credentials file.  Default: ~/opt/config/gdocs_credentials")

    parser.add_option("-c", "--conf", dest="conf", 
    default=os.path.join(os.environ['HOME'],'opt/config/post_process.yaml'), 
    help = "Config file.  Default: ~/opt/config/post_process.yaml")

    (options, args) = parser.parse_args()

    LOG = scilifelab.log.file_logger('LOG', options.conf ,'analysis2coucdb.log')
    main(options.project_name, options.conf, options.cred)

