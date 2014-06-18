#!/home/hiseq.bioinfo/.virtualenvs/lims2db/bin/python
from __future__ import print_function
"""Script to load project info from Lims into the project database in statusdb.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""

import sys
import os
import codecs
from optparse import OptionParser
from scilifelab.db.statusDB_utils import *
from helpers import *
from pprint import pprint
from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD
import objectsDB as DB
from datetime import date
import scilifelab.log
lims = Lims(BASEURI, USERNAME, PASSWORD)

def print_couchdb_obj_to_file(obj, output_f=None):
    if output_f is not None:
        with open(output_f, 'w') as f:
            print(obj, file=f)
    else:
        print(obj, file=sys.stdout)

def  main(proj_name, all_projects, days, conf, upload_data=True, output_f=None):
    first_of_july = '2013-06-30'
    today = date.today()
    couch = load_couch_server(conf)
    proj_db = couch['projects']
    samp_db = couch['samples']
    if all_projects:
        projects = lims.get_projects()
        for proj in projects:
            LOG.info(proj.name)
            try:
                closed = date(*map(int, proj.close_date.split('-')))
                days_closed = (today-closed).days
            except:
                days_closed = 0
            if proj.open_date:
                ordered_opened = proj.open_date
            elif 'Order received' in dict(proj.udf.items()).keys():
                ordered_opened = proj.udf['Order received'].isoformat()
            else:
                LOG.info("Project is not updated because 'Order received' date"
                    " and 'open date' is missing for project %s" % proj.name)
            if ordered_opened:
                if comp_dates(first_of_july, ordered_opened):
                    if (days_closed < days):
                        proj_time = time.time()
                        try:
                            obj = DB.ProjectDB(lims, proj.id, samp_db)
                            key = find_proj_from_view(proj_db, proj.name)
                            obj.project['_id'] = find_or_make_key(key)
                            if upload_data:
                                info = save_couchdb_obj(proj_db, obj.project)
                            else:
                                info = print_couchdb_obj_to_file(obj.project, 
                                    output_f)
                            LOG.info('project %s is handeled and %s : _id = %s' 
                                % (proj.name, info, obj.project['_id']))
                        except:
                            LOG.info('Issues geting info for %s. The '
                                '"Application" udf might be missing' % proj.name)
                    else:
                        LOG.info('Project is not updated because the project '
                            'has been closed for %s days.' % days_closed)  
                else:
                    LOG.info(('Project is not updated because the project was'
                        ' opened or ordered before 2013-06-30 (%s)') % ordered_opened)
    elif proj_name is not None:
        proj = lims.get_projects(name = proj_name)
        if len(proj) == 0:
            LOG.warning('No project named %s in Lims' % proj_name)
        else:
            proj = proj[0]
            if proj.open_date:
                ordered_opened = proj.open_date
            elif 'Order received' in dict(proj.udf.items()).keys():
                ordered_opened = proj.udf['Order received'].isoformat()
            else:
                LOG.info("Project is not updated because 'Order received' date"
                    " and 'open date' is missing for project %s" % proj.name)
            if ordered_opened:
                if comp_dates(first_of_july, ordered_opened):
                    cont = 'yes'
                else:
                    cont = raw_input('The project %s is ordered before '
                        '2013-07-01. Do you still want to load the data from '
                        'lims into statusdb? (yes/no): ' % proj_name)
                if cont == 'yes':
                    obj = DB.ProjectDB(lims, proj.id, samp_db)
                    key = find_proj_from_view(proj_db, proj.name)
                    obj.project['_id'] = find_or_make_key(key)
                    if upload_data:
                        info = save_couchdb_obj(proj_db, obj.project)
                    else:
                        info = print_couchdb_obj_to_file(obj.project, output_f)
                    LOG.info('project %s is handeled and %s : _id = %s' 
                                        % (proj_name, info, obj.project['_id']))

if __name__ == '__main__':
    usage = "Usage:       python project_summary_upload_LIMS.py [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-p", "--project", dest = "project_name", default = None,
        help = "eg: M.Uhlen_13_01. Dont use with -a flagg.")

    parser.add_option("-a", "--all_projects", dest = "all_projects",
        action = "store_true", default = False, help = ("Upload all Lims ",
        "projects into couchDB. Don't use with -f flagg."))

    parser.add_option("-d", "--days", dest = "days", default = 60, help = (
        "Projects with a close_date older than DAYS days are not updated. ",
        "Default is 60 days. Use with -a flagg"))

    parser.add_option("-c", "--conf", dest = "conf", default = os.path.join(
        os.environ['HOME'],'opt/config/post_process.yaml'), help = ("Config ",
        "file.  Default: ~/opt/config/post_process.yaml"))

    parser.add_option("--no_upload", dest = "upload", default = True, 
        action = "store_false", help = ("Use this tag if project objects ",
        "should not be uploaded, but printed to output_f, or to stdout"))

    parser.add_option("--output_f", dest = "output_f", help = ("Output file",
        " that will be used only if --no_upload tag is used"))

    (options, args) = parser.parse_args()
    LOG = scilifelab.log.file_logger('LOG',options.conf,
        'lims2db_projects-tools.log', 'log_dir_tools')
    main(options.project_name, options.all_projects, options.days, options.conf,
         upload_data = options.upload, output_f = options.output_f)

