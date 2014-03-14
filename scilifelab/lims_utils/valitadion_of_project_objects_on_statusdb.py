#!/usr/bin/env python

"""valitadion_of_LIMS_upgrade.py is a script to compare extraction output from lims stage 
server and lims production server. The comparison is based on the objects created to build 
documents in the projects database on status db. A recursive function compares all values 
in the objects and any differing values or missing keys are logged in a validation log file.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""

usage = """

*****Recomended validation procedure:*****

Testing the script:

Test that the script is caching differences by changing something on the 
stage server, eg. the value of the sample udf "status_(manual)". for some 
project J.Doe_00_00. Then run the script with the -p flagg:

valitadion_of_LIMS_upgrade.py -p J.Doe_00_00

This should give the output:

Lims stage and Lims production are differing for proj J.Doe_00_00: True
Key status_(manual) differing: Lims production gives: Aborted. Lims stage gives In Progress.

Running the validation:

Run valitadion_of_LIMS_upgrade.py with the -a flagg and grep for "True" in 
the logfile when the script is finished. It will take some hours to go through 
all projects opened after jul 1

If you don't find anything when grepping for True in the log file, no differences 
are found for any projects.

If you get output when grepping for True, there are differences. Then read the log 
file to find what is differing. 

"""
import sys
import os
import codecs
from optparse import OptionParser
from scilifelab.db.statusDB_utils import *
from helpers import *
from pprint import pprint
import objectsDB as DB
from datetime import date
import scilifelab.log

def comp_obj(stage, prod):
    """compares the two dictionaries obj and dbobj"""
    LOG.info('project %s is being handeled' % stage['project_name'])
    diff = recursive_comp(stage, prod)
    LOG.info('tools and tools-dev are differing for proj %s: %s' % ( stage['project_name'],diff))

def recursive_comp(stage, prod):
    diff = False
    keys = list(set(stage.keys() + prod.keys()))
    for key in keys:
        if not (stage.has_key(key)):
            LOG.info('Key %s missing in tools to db object ' % key)
            diff = True
        elif not  prod.has_key(key):
            LOG.info('Key %s missing in tools-dev to db object ' % key)
            diff = True
        else:
            prod_val = prod[key]
            stage_val = stage[key]
            if (prod_val != stage_val):
                diff = True
                if (type(prod_val) is dict) and (type(stage_val) is dict):
                    diff = diff and recursive_comp(stage_val, prod_val)
                else:
                    LOG.info('Key %s differing: tools gives: %s. tools-dev gives %s. ' %( key,prod_val,stage_val))
    return diff

def  main(proj_name, all_projects, conf_tools_dev):
    couch_tools = load_couch_server(os.path.join(os.environ['HOME'],'opt/config/post_process.yaml'))
    couch_tools_dev = load_couch_server(conf_tools_dev)
    proj_db_tools = couch_tools['projects']
    proj_db_tools_dev = couch_tools_dev['projects']
    if all_projects:
        for key in proj_db_tools:
            proj_tools = proj_db_tools.get(key)
            proj_tools_dev = proj_db_tools_dev.get(key)
            proj_name = proj_tools['project_name']
            try:
                if not proj_tools_dev:
                    LOG.warning("""Found no projects on tools-dev with name %s""" % proj_name)
                else:
                    comp_obj(proj_tools, proj_tools_dev)
            except:
                LOG.info('Failed comparing stage and prod for proj %s' % proj_name)    
    elif proj_name is not None:
        key = find_proj_from_view(proj_db_tools, proj_name)
        proj_tools = proj_db_tools.get(key)
        proj_tools_dev = proj_db_tools_dev.get(key)
        if (not proj_tools) | (not proj_tools_dev):
            LOG.warning("Found no project named %s" %(proj_name))
        else:
            comp_obj(proj_tools, proj_tools_dev)

if __name__ == '__main__':
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--project", dest="project_name", default=None,
    help = "eg: M.Uhlen_13_01. Dont use with -a flagg.")

    parser.add_option("-a", "--all_projects", dest="all_projects", action="store_true", default=False,
    help = "Check all projects on couchDB. Don't use with -f flagg.")

    parser.add_option("-c", "--conf", dest="conf", 
    default=os.path.join(os.environ['HOME'],'opt/scilifelab/scilifelab/lims_utils/post_process.yaml'),         
    help = "Config file for tools-dev.  Default: ~/opt/scilifelab/scilifelab/lims_utils/post_process.yaml")

    (options, args) = parser.parse_args()

    LOG = scilifelab.log.file_logger('LOG',options.conf ,'validate_projects_statusdb.log', 'log_dir_tools')
    main(options.project_name, options.all_projects, options.conf)

