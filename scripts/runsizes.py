#!/usr/bin/env python
'''
Gets filesize of a first level of directories and sends it to a CouchDB instance.

Not using os.path.getsize neither os.stat.st_size since they report
inaccurate filesizes:

http://stackoverflow.com/questions/1392413/calculating-a-directory-size-using-python

Usage: python runsizes.py --server 'http://localhost:5984' --db size_logs_tests --dir .

'''
# TODO: Manage depth of root
# TODO: Filter out by .bcl files and/or include other ones

import os
import sys
import argparse
import subprocess
import datetime
import couchdb


def get_dirsizes(path="."):
    ''' Gets directory size.
        TODO: Be replaced with a more pythonic way which reports the size correctly.
    '''
    path = path.strip().split('\t')
    out = subprocess.check_output(["du", "-s", path[0]])
    return out.split('\t')[0]


def main():
    dirsizes = {"time": datetime.datetime.now().isoformat()}
    server = 'localhost:5984'
    db = "tests"
    root = "."

    parser = argparse.ArgumentParser(description="Compute directory size(s) and report them to a CouchDB database")

    parser.add_argument('--dir', dest='root', action='store',
                        help="the directory to calculate dirsizes from")

    parser.add_argument("--server", dest='server', action='store',
                        help="CouchDB instance to connect to, defaults to localhost:5984")

    parser.add_argument("--db", dest='db', action='store',
                       help="CouchDB database name, defaults to 'tests'")

    args = parser.parse_args()

    for d in os.listdir(args.root):
    	path = os.path.join(args.root, d)
        dirsizes[path] = get_dirsizes(path)

    couch = couchdb.Server(args.server)
    db = couch[args.db]
    db.save(dirsizes)

if __name__ == "__main__":
    main()
