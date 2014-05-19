#! /usr/bin/env python

import argparse
import ConfigParser
import datetime
import getpass
import os
import platform

import couchdb
from crontab import CronTab
from couchdb import ResourceNotFound

def update_cronjobs_database(couch):
    """ Update cronjobs database with information of the users' crontab

    :param couch: CouchDB Server instance
    """
    crontab = CronTab()
    server = platform.node().split('.')[0]
    try:
        cdb = couch['cronjobs']
        server_view = cdb.view('server/alias')
        crontab_json = {'server': server, 'cronjobs': {}, 'Last updated': str(datetime.datetime.now())}

        for job_id, job in enumerate(crontab.crons):
            job_json = {}
            job_json['User'] = getpass.getuser()
            job_json['Command'] = job.command
            job_json['Enabled'] = job.enabled
            job_json['Comment'] = job.comment
            job_json['Minute'] = str(job.minute)
            job_json['Hour'] = str(job.hour)
            job_json['Day of month'] = str(job.dom)
            job_json['Month'] = str(job.month)
            job_json['Day of week'] = str(job.day)
            crontab_json['cronjobs'][job_id] = job_json
                 
        # There can only be one document per server, if so
        server_doc = {}
        for row in server_view.rows:
            if row.key == server:
                server_doc = cdb.get(row.value)
        if server_doc:
            server_doc.update(crontab_json)
        else:
            server_doc = crontab_json
        cdb.save(server_doc)
        
    except ResourceNotFound as e:
        raise e('ERROR while looking for cronjobs database in {}'.format(couch.resource.url))

if __name__=='__main__':
    DESCRIPTION = """ Script to update StatusDB with crontab information.

    Will pick up information from the crontab and update it to StatusDB. The intention
    of this is that everyone can track what tasks are running on the background and
    when.

    Requires module python-crontab: pip install python-crontab
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--config', type=str, help="Configuration file with StatusDB " + 
        "credentials. By default will try to load ~/.pm/pm.conf")
    args = parser.parse_args()

    config = ConfigParser.SafeConfigParser()
    config_file = args.config if args.config else os.path.join(os.environ['HOME'], '.pm/pm.conf')
    try:
        with open(config_file, 'r') as f:
            config.readfp(f)

        db = config.get('db', 'url').rstrip()
        user = config.get('db', 'user').rstrip()
        password = config.get('db', 'password').rstrip()
        port = '5984' if not config.has_option('db', 'port') else config.get('db', 'port')
        couch = couchdb.Server('http://' + ':'.join([db, port]))
        couch.resource.credentials = (user, password)
        update_cronjobs_database(couch)
    except IOError as e:
        raise e("Please make sure you've created your own configuration file " +  
            "(i.e: ~/.pm/pm.conf), and that it contains a source and a destination servers")
