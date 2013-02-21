#!/usr/bin/env python
"""Gets filesize of a first level of directories and sends it to a CouchDB instance.

Not using os.path.getsize neither os.stat.st_size since they report
inaccurate filesizes:

http://stackoverflow.com/questions/1392413/calculating-a-directory-size-using-python
"""
# TODO: Manage depth of root (how many dir levels): http://stackoverflow.com/questions/229186/os-walk-without-digging-into-directories-below
# TODO: Filter out by .bcl files and/or include other ones

import os
import argparse
import subprocess
import datetime
import couchdb
import re


def get_dirsizes(path="."):
    """Gets directory size.
        TODO: Be replaced with a more pythonic way which reports the size correctly.
    """
    path = path.strip().split('\t')
    out = subprocess.check_output(["du", "-sb", path[0]], stderr=subprocess.STDOUT)
    return out.split('\t')[0]


def parse_dirsizes(path, dirsizes={"errors": []}):
    """Parse directory sizes that have been saved to a file
    """

    date_regexp = r'(?:\d{2})?(?:\d{2}\-?){3}[\sT]\d{2}(?:\:\d{2}){1,2}'
    try:
        with open(path) as fh:
            for line in fh:
                # Parse a timestamp
                m = re.search(date_regexp, line)
                if m:
                    try:
                        timestamp = datetime.datetime.strptime(m.group(0), "%Y-%m-%d %H:%M")
                        dirsizes["time"] = timestamp.isoformat()
                    except:
                        pass

                    continue

                # Assume directories are listed as [size] [path] on one line each
                try:
                    splits = line.split()
                    if len(splits) < 2:
                        continue

                    size, path = int(splits[0]), splits[1]
                    dirsizes[path] = size
                except ValueError:
                    continue

    except Exception as e:
        dirsizes["errors"].append(str(e))

    return dirsizes


def send_db(server, db, data):
    """ Submits provided data to database on server
    """
    couch = couchdb.Server(server)
    db = couch[db]
    db.save(data)
    #with open("runsizes.log", "w") as fh:
    #   print "Saving data to %s" % fh
    #   fh.write(str(_to_unicode(data)))
    #   print "Sending data to couchdb"


def main():
    dirsizes = {"time": datetime.datetime.now().isoformat(),
                "unit": "bytes",
                "errors": []}

    parser = argparse.ArgumentParser(description="Compute directory size(s) and report them to a CouchDB database")

    parser.add_argument('--dir', dest='root', action='append',
                        help="the directory to calculate dirsizes from")

    parser.add_argument("--server", dest='server', action='store', default="localhost:5984",
                        help="CouchDB instance to connect to, defaults to localhost:5984")

    parser.add_argument("--db", dest='db', action='store', default="tests",
                        help="CouchDB database name, defaults to 'tests'")

    parser.add_argument("--dry-run", dest='dry_run', action='store_true', default=False,
                        help="Do not submit the resulting hash to CouchDB")

    args = parser.parse_args()

    for r in args.root:  # multiple --dir args provided
        if os.path.exists(r) and os.path.isdir(r):
            for d in os.listdir(r):
                path = os.path.join(r, d)
                try:
                    dirsizes[path] = int(get_dirsizes(path))
                except subprocess.CalledProcessError as pe:
                    dirsizes['errors'].append(pe.output)
        else:
            dirsizes = parse_dirsizes(r, dirsizes)

    if args.dry_run:
        print(dirsizes)
    else:
        send_db(args.server, args.db, dirsizes)


if __name__ == "__main__":
    main()
