#!/usr/bin/env python
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
import objectsDB as DB_v0
from datetime import date
import time
import scilifelab.log
lims = Lims(BASEURI, USERNAME, PASSWORD)

class PSUL():
    def __init__(self, proj, samp_db, proj_db, upload_data, days, man_name):
        self.proj = proj
        self.id = proj.id
        self.udfs = proj.udf
        self.name = proj.name
        self.open_date = proj.open_date
        self.close_date = proj.close_date
        self.samp_db = samp_db
        self.proj_db = proj_db
        self.upload_data = upload_data
        self.man_name = man_name
        self.days = days

    def print_couchdb_obj_to_file(self, obj, output_f = None):
        if output_f is not None:
            with open(output_f, 'w') as f:
                print(obj, file = f)
        else:
            print(obj, file = sys.stdout)

    def get_ordered_opened(self):
        """Is project registered as opened or ordered?"""

        if self.open_date:
            return self.open_date
        elif 'Order received' in dict(self.udfs.items()).keys():
            return self.udfs['Order received'].isoformat()
        else:
            LOG.info("Project is not updated because 'Order received' date and "
                     "'open date' is missing for project %s" % self.name)
            return None

    def get_days_closed(self):
        """Project registered as closed?"""

        if self.close_date:
            closed = date(*map(int, self.close_date.split('-')))
            return (date.today() - closed).days
        else:
            return 0

    def determine_update(self,ordered_opened):
        """Determine wether to and how to update project"""
        days_closed = self.get_days_closed()
        opende_after_130630 = comp_dates('2013-06-30', ordered_opened)
        opende_after_140630 = comp_dates('2014-06-30', ordered_opened)
        closed_for_a_while = (days_closed > self.days)
        log_info = ''
        LOG.info(days_closed)
        LOG.info(self.days)
        LOG.info(closed_for_a_while)
        LOG.info(ordered_opened)
        if (not opende_after_130630) or closed_for_a_while:
            if self.man_name:   ## Ask wether to update
                start_update = raw_input("""
                Project %s was ordered or opended at %s and has been closed for 
                %s days. Do you still want to load the data from lims into statusdb? 
                Press enter for No, any other key for Yes!  """ % (self.name, 
                ordered_opened, days_closed))
            else:               ## Do not update
                start_update = False
                log_info = ('Project is not updated because: ')
                if closed_for_a_while:
                    log_info += ('It has been closed for %s days. ' % days_closed)
                if not opende_after_130630:
                    log_info += ('It was opened or ordered before 2013-06-30 (%s)' 
                                 % ordered_opened)
        else:
            start_update = True

        if start_update:
            database = DB_v1 if opende_after_140630 else DB_v0
            return log_info, database
        else:
            return log_info, None


    def update_project(self, database):
        """Fetch project info and update project in the database."""
        try:
            obj = database.ProjectDB(lims, self.id, self.samp_db)
            key = find_proj_from_view(self.proj_db, self.name)
            obj.project['_id'] = find_or_make_key(key)
            if self.upload_data:
                info = save_couchdb_obj(self.proj_db, obj.project)
            else:
                info = self.print_couchdb_obj_to_file(obj.project, output_f)
            return ('project %s is handeled and %s : _id = %s' %(self.name,
                            info, obj.project['_id']))
        except:
            return ('Issues geting info for %s. The "Application" udf '
                                'might be missing' % self.name)

    def project_update_and_logging(self, proj_num = '', num_projs = ''):
        start_time = time.time()
        ordered_opened = self.get_ordered_opened()
    
        if ordered_opened:
            log_info, database = self.determine_update(ordered_opened)
            if database:
                log_info = self.update_project(database)
        else:
            log_info = ('No open date or order date found for project %s. '
                        'Project not updated.' % self.name) 

        elapsed = time.time() - start_time
        LOG.info('Time - %s : Proj - %s/%s : Name - %s' % (elapsed, proj_num, 
                                                        num_projs , self.name))
        LOG.info(log_info) 

def main(man_name, all_projects, days, conf, upload_data, output_f = None):
    couch = load_couch_server(conf)
    proj_db = couch['projects']
    samp_db = couch['samples']

    if all_projects:
        projects = lims.get_projects()
        num_projs = len(projects)
        for proj_num, proj in enumerate(projects):
            P = PSUL(proj, samp_db, proj_db, upload_data, days, man_name)
            P.project_update_and_logging(proj_num, num_projs)
    elif man_name:
        proj = lims.get_projects(name = man_name)
        if not proj:
            LOG.warning('No project named %s in Lims' % man_name)
        else:
            P = PSUL(proj[0], samp_db, proj_db, upload_data, days, man_name)
            P.project_update_and_logging()


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
    LOG = scilifelab.log.file_logger('LOG', options.conf, 'lims2db_projects.log'
                                                               ,'log_dir_tools')
    main(options.project_name, options.all_projects, options.days, options.conf,
         upload_data = options.upload, output_f = options.output_f)

