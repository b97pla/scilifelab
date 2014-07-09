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

def comp_obj(proj_tools_dev, diff, keys):
    """compares the two dictionaries obj and dbobj"""
    if proj_tools_dev.has_key('project_name'):
        if keys:
            print keys
            LOG.info('tools and tools-dev are differing for proj %s: %s' % (
                                                proj_tools_dev['project_name'],diff))
        LOG.info('tools and tools-dev are differing for proj %s: %s' % ( 
                                    proj_tools_dev['project_name'],diff))
    elif proj_tools_dev:
        LOG.info('project_name missing in %s' % (proj_tools_dev['_id']))

def missing_keys(proj_tools_dev, proj_tools,proj_name,application,diff_keys=[],diff = False):
    G20158 = ['m_reads_sequenced','status','size_(bp)']
    BASEDIF = ['_rev','modification_time','reagent_labels']
    keys = list(set(proj_tools_dev.keys() + proj_tools.keys()))
    for key in keys:
        if len(key.split('_'))==2 and key[0]=='P':
                        LOG.info('SAMPLE:   %s' %(key))
        if (key not in BASEDIF+G20158) and (proj_tools.has_key(key)) and proj_tools[key]:
            if proj_tools_dev.has_key(key) and proj_tools_dev[key]:
                proj_tools_val = proj_tools[key]
                proj_tools_dev_val = proj_tools_dev[key]
                if (proj_tools_val != proj_tools_dev_val):
                    diff = True
                    if (type(proj_tools_val) is dict) and (
                                        type(proj_tools_dev_val) is dict):
                        diff,diff_keys = missing_keys(proj_tools_dev_val, 
                               proj_tools_val,proj_name,application,diff_keys=diff_keys,diff = diff)
            elif not proj_tools_dev.has_key(key):
                LOG.info('tools-dev proj %s-%s - missing Key %s' %(proj_name,
                                                application,key))
                diff = True
                diff_keys.append(key)
    return diff, diff_keys


def differing_keys(proj_tools_dev,proj_tools,proj_name,application,diff_keys = [],diff = False):
    G20158 = ['m_reads_sequenced','status','size_(bp)']
    BASEDIF = ['_rev','modification_time']
    keys = list(set(proj_tools_dev.keys() + proj_tools.keys()))
    for key in keys:
        if len(key.split('_'))==2 and key[0]=='P':
            LOG.info('SAMPLE:   %s' %(key))
        if key not in BASEDIF+G20158:
            if proj_tools.has_key(key) and proj_tools_dev.has_key(key):
                if not proj_tools_dev[key]:
                    LOG.info('empty key %s in proj %s'%(key,proj_name))
                proj_tools_val = proj_tools[key]
                proj_tools_dev_val = proj_tools_dev[key]
                if (proj_tools_val != proj_tools_dev_val):
                    if (type(proj_tools_val) is dict) and (type(
                                                   proj_tools_dev_val) is dict):
                        diff, diff_keys = differing_keys(proj_tools_dev_val, 
                               proj_tools_val,proj_name,application, diff_keys,diff = diff)
                    else:
                        diff = True
                        diff_keys.append(key)
                        LOG.info('Key %s is differing for proj %s-%s: tools gives:' 
                                 '%s. tools-dev gives %s. ' %( key, proj_name,
                                 application,proj_tools_val,
                                 proj_tools_dev_val))
    return diff, diff_keys

def  main(proj_name, all_projects, conf_tools_dev, conf_tools):
    couch_tools = load_couch_server(conf_tools)
    couch_tools_dev = load_couch_server(conf_tools_dev)
    proj_db_tools = couch_tools['projects']
    proj_db_tools_dev = couch_tools_dev['projects']
    if all_projects:
        proj_name = ''
        application = ''
        for key in proj_db_tools:
            keys =[]
            proj_tools = proj_db_tools.get(key)
            proj_tools_dev = proj_db_tools_dev.get(key)
            try:
                proj_name = proj_tools['project_name']
                LOG.info('Handeling %s %s' % (proj_name, key))
                try:
                    application = proj_tools['application']
                except:
                    LOG.info('Application is missing for proj %s in tools' %(proj_name))
            except:
                proj_name = 'XX'
                LOG.info('Project name missing in tools: %s' %(key))
            if not proj_tools_dev:
                LOG.warning("Found no projects on tools-dev with name %s" % 
                                                                  proj_name)
            elif proj_name and application:
                diff , keys= differing_keys(proj_tools_dev, proj_tools, proj_name, 
                                                                application, diff_keys= keys)
                print diff
                diff, keys = missing_keys(proj_tools_dev, proj_tools, proj_name,
                                                    application,diff_keys=keys,diff = diff)
                print diff
                comp_obj(proj_tools_dev, diff,keys)
    elif proj_name is not None:
        keys =[]
        key = find_proj_from_view(proj_db_tools, proj_name)
        proj_tools = proj_db_tools.get(key)
        proj_tools_dev = proj_db_tools_dev.get(key)
        LOG.info('Handeling %s %s' % (proj_name, key))
        if (not proj_tools) | (not proj_tools_dev):
            LOG.warning("Found no project named %s" %(proj_name))
        else:
            diff,keys = differing_keys(proj_tools_dev, proj_tools,proj_name,
                                             proj_tools['application'], diff_keys=keys)
            print diff
            diff,keys = missing_keys(proj_tools_dev, proj_tools, proj_name,
                                         proj_tools['application'],diff_keys=keys,diff = diff)
            print diff
            comp_obj(proj_tools_dev, diff,keys)

if __name__ == '__main__':
    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--project", dest="project_name", default=None,
        help = "eg: M.Uhlen_13_01. Dont use with -a flagg.")
    parser.add_option("-a", "--all_projects", dest="all_projects", 
        action="store_true", default=False,
        help = "Check all projects on couchDB. Don't use with -f flagg.")
    parser.add_option("-c", "--conf", dest="conf", 
        help = "Config file for tools.")
    parser.add_option("-d", "--conf_dev", dest="conf_dev", 
        help = "Config file for tools-dev.") 
    (options, args) = parser.parse_args()

    LOG = scilifelab.log.file_logger('LOG',options.conf ,
                                     'validate_projects_statusdb.log', 
                                     'log_dir_tools')
    main(options.project_name, options.all_projects, options.conf_dev, 
                                                        options.conf)

