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
LOG = scilifelab.log.minimal_logger('LOG')

def  main(flowcell, all_flowcells,days,conf):
    """If all_flowcells: all runs run less than a moth ago are uppdated"""
    today = date.today()
    couch = load_couch_server(conf)
    fc_db = couch['flowcells']
    if all_flowcells:
        flowcells = lims.get_processes(type = ['Illumina Sequencing (Illumina SBS) 4.0','MiSeq Run (MiSeq) 4.0'])
        for fc in flowcells:
            try:
                closed = date(*map(int, fc.date_run.split('-')))
                delta = today-closed
            except AttributeError:
                #Happens if fc has no date run, we should just not update and get to the next flowcell
                continue

            #if delta.days < days and dict(fc.udf.items()).has_key('Flow Cell ID'):
            if dict(fc.udf.items()).has_key('Flow Cell ID'):
                if '-' in dict(fc.udf.items())['Flow Cell ID']:
                    flowcell_name = dict(fc.udf.items())['Flow Cell ID']
                elif dict(fc.udf.items()).has_key('Flow Cell Position'):
                    flowcell_name = dict(fc.udf.items())['Flow Cell Position'] + dict(fc.udf.items())['Flow Cell ID'] 
                key = find_flowcell_from_view(fc_db, flowcell_name)
                if key:
                    dbobj = fc_db.get(key)
                    print dbobj['modification_time']+'  '+key
                    if delta.days < days:
                        dbobj["illumina"]["run_summary"] = get_sequencing_info(fc)
                        info = save_couchdb_obj(fc_db, dbobj)
                        LOG.info('flowcell %s %s : _id = %s' % (flowcell_name, info, key))
    elif flowcell is not None:
        if '-' in flowcell:
            flowcell_name = flowcell
            fc = lims.get_processes(type = 'MiSeq Run (MiSeq) 4.0', udf = {'Flow Cell ID' : flowcell_name})[0]
        else:
            flowcell_name_short = flowcell[1:len(flowcell)]
            flowcell_position = flowcell[0]
            fc = lims.get_processes(type = 'Illumina Sequencing (Illumina SBS) 4.0', udf = {'Flow Cell ID' : flowcell_name_short, 'Flow Cell Position' : flowcell_position})[0]
            flowcell_name = flowcell_position + flowcell_name_short
        key = find_flowcell_from_view(fc_db, flowcell_name)
        if key:
            dbobj = fc_db.get(key)
            dbobj["illumina"]["run_summary"] = get_sequencing_info(fc)
            get_run_qcs(fc, dbobj['lanes'])
            info = save_couchdb_obj(fc_db, dbobj)
            LOG.info('flowcell %s %s : _id = %s' % (flowcell_name, info, key))
                

if __name__ == '__main__':
    usage = "Usage:       python flowcell_summary_upload_LIMS.py [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-f", "--flowcell", dest="flowcell_name", default=None, 
    help = "eg: AD1TAPACXX. Don't use with -a flagg.")

    parser.add_option("-a", "--all_flowcells", dest="all_flowcells", action="store_true", default=False, 
    help = "Uploads all Lims flowcells into couchDB. Don't use with -f flagg.")

    parser.add_option("-d", "--days", dest="days", default=30, 
    help="Runs older than DAYS days are not updated. Default is 30 days. Use with -a flagg")

    parser.add_option("-c", "--conf", dest="conf", 
    default=os.path.join(os.environ['HOME'],'opt/config/post_process.yaml'), 
    help = "Config file.  Default: ~/opt/config/post_process.yaml")

    (options, args) = parser.parse_args()

    LOG = scilifelab.log.file_logger('LOG', options.conf, 'lims2db_flowcells.log','log_dir_tools')
    main(options.flowcell_name, options.all_flowcells, options.days, options.conf)

