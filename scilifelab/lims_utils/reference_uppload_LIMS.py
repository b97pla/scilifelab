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
from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD
from datetime import date
from lims_utils import *
from scilifelab.db.statusDB_utils import *
import scilifelab.log
lims = Lims(BASEURI, USERNAME, PASSWORD)

def  main(project, sample, conf):
    today = date.today()
    couch = load_couch_server(conf)
    ref_db = couch['reference']
    if project:
        load_project_udfs(ref_db)
    if sample:
        load_sample_udfs(ref_db)

def load_sample_udfs(ref_db):
    """Loads all sample level udfs into the reference database on statusdb.
    ref_db is the reference database."""
    udfs = prepare_udf_source_info('Sample', lims.get_samples())    
    samples = ref_db.get('project-samples-[key]')
    samples_details = ref_db.get('project-samples-[key]-details')
    for udf, udf_field in udfs.items():
        if 'discontinued' not in udf.split('_'):
            if udf in SAMP_UDF_EXCEPTIONS:
                key = 'project-samples-[key]-'+udf
                samples[udf] = key
                inf=save_couchdb_ref_obj(ref_db, samples)
            else:
                key = 'project-samples-[key]-details-'+udf
                samples_details[udf] = key
                delkey = 'project-samples-[sample]-details-'+udf
                delete_doc(ref_db,delkey)
                inf=save_couchdb_ref_obj(ref_db, samples_details)
            udf_field['_id'] = key
            inf=save_couchdb_ref_obj(ref_db, udf_field)

def load_project_udfs(ref_db):
    """Loads all project level udfs into the reference database on statusdb.
    ref_db is the reference database."""
    udfs = prepare_udf_source_info('Project', lims.get_projects())
    project = ref_db.get('project')
    project_details = ref_db.get('project-details')
    for udf, udf_field in udfs.items():
        if 'discontinued' not in udf.split('_'):
            if udf in PROJ_UDF_EXCEPTIONS:
                key = 'project-'+udf
                project[udf] = key
                inf=save_couchdb_ref_obj(ref_db, project)
            else:
                key = 'project-details-'+udf
                project_details[udf] = key
                inf=save_couchdb_ref_obj(ref_db, project_details)
            udf_field['_id'] = key
            inf=save_couchdb_ref_obj(ref_db, udf_field)

def prepare_udf_source_info(element_type = None, element_list = None):
    """element_type should be some lims element type, such as 
    'Project', 'Sample', 'Artifact'... 
    from wich you want to get udf source info.
    element_list should be an list of instanses of the element type. 
    eg output from get_samples(), get_projects()..."""
    udfs = lims.get_udfs(attach_to_name = element_type)
    objects={}
    for udf in udfs:
        db_name=udf.name.lower().replace(' ','_').replace('-','_')
        objects[db_name] = {'doc_source': { 'lims_field':udf.name,
                                            'lims_element': element_type, 
                                            'source': 'Lims'},
                            'doc_type': udf.root.get('type')}
    return objects  

def find_example(element_list = None, udf_name = None):
    """To find an example value od the udf_name in the element_list, 
    where element_list should be an list of instanses of some element type. 
    Eg. output from get_samples(), get_projects()... Will only make 1000 tryals."""
    i=0
    if element_list and udf_name:
        print 'Searching udf-example of: '+ udf_name
        for elm in element_list:
            i=i+1
            try:
                example=dict(elm.udf.items())[udf_name]
                if example:
                    print example
                    return example
            except:
                pass 
            if i==1000:
                return 'Not found'

if __name__ == '__main__':
    usage = "Usage:       python flowcell_summary_upload_LIMS.py [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--project", dest="project", action="store_true", default=False, 
    help = "Upload source info for project level udfs into the reference database.")
    
    parser.add_option("-s", "--sample", dest="sample", action="store_true", default=False,            
    help = "Upload source info for sample level udfs into the reference database.")

    parser.add_option("-c", "--conf", dest="conf", 
    default=os.path.join(os.environ['HOME'],'opt/config/post_process.yaml'), 
    help = "Config file.  Default: ~/opt/config/post_process.yaml")

    (options, args) = parser.parse_args()

    LOG = scilifelab.log.file_logger('LOG', options.conf, 'lims2db_reference.log')
    main(options.project, options.sample , options.conf)

