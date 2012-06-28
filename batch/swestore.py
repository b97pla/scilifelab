#!/usr/bin/env python
'''
Stores datasets that are older than --period into long term storage and reports
back to statusDB
'''

import os
import sys
import argparse
import subprocess
import datetime
import couchdb
import drmaa
import StringIO
import logging
import multiprocessing

def send_db(server, db, data):
    ''' Submits provided data to database on server
    '''
    couch = couchdb.Server(server)
    db = couch[db]
    db.save(data)

def submit_drmaa(cmd, time, part, prj, name):
    ''' Submits job to compress dataset
    '''
    jobid = None

    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = cmd[0]
    jt.args = cmd[1:]
    jt.nativeSpecification = "-A %s -p %s -t %s -J %s" % (prj, part, time, name)
    
    try:
        jobid = s.runJob(jt)
        status = s.jobStatus(jobid)
    except drmaa.errors.InternalException as e:
        s.deleteJobTemplate(jt)
        s.exit()
        print e.message

    return jobid, status

def find_older_than(path, days):
    ''' Wrapper function for find functionality
        WARNING: Tailored for illumina specific flowcell glob: "*_*_*_*"

        XXX: Find a more pythonic way for this
    ''' 
    out = None
    dsets = []

    try:
        out = subprocess.check_output(["find", path[0], "-maxdepth", "1", "-iname", \
                                       "*_*_*_*", "-type", "d", "-mtime", '+'+str(days)], \
                                        stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as pe:
        print pe.output

    out = out.rstrip()
    dsets = out.split("\n")

    return dsets


def main():
    logstamp = {"time": datetime.datetime.now().isoformat(),
                "jobid": None,
                "status": None
    }
    
    datasets = {}


    parser = argparse.ArgumentParser(description="Archive datasets and report status")

    parser.add_argument('--src-dir', dest='root', action='append',
                        help="The directory containing datasets")

    parser.add_argument('--dst-dir', dest='dst', action='store',
                        help="Destination directory (to be picked and archived to long-term storage)")
    
    parser.add_argument('--days', dest='days', action='store',
                        help="Number of days to keep the datasets on primary storage (away from longterm)")

    parser.add_argument('--server', dest='server', action='store', default='localhost:5984',
                        help="CouchDB instance to connect to, defaults to localhost:5984")

    parser.add_argument('--db', dest='db', action='store', default='tests',
                        help="CouchDB database name, defaults to 'tests'")

    parser.add_argument('--dry-run', action='store_true', default=False,
                        help="Do not submit the resulting hash to CouchDB")

    args = parser.parse_args()


    # XXX check DRMAA status, report to db only if it changes.
    # decouple that into another script ?

    datasets = find_older_than(args.root, args.days)

    if not args.dry_run:
        for dset in datasets:
            archive_to = args.dst+os.path.basename(dset)+str(".tar.bz2")
            cmd = ["tar", "cf", archive_to, '--use-compress-prog=pbzip2', dset]
            
            jobid, status = submit_drmaa(cmd, "00:00:10", "node", "a2010002", "swestore")
            logstamp["jobid"] = int(jobid)
            logstamp["status"] = status
       
            #send_db(args.server, args.db, logstamp)

if __name__ == "__main__":
    main()
