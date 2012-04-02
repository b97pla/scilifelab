#!/usr/bin/env python
'''
Gets filesize of a first level of directories and reports
it in JSON format.
Not using os.path.getsize neither os.stat.st_size since they report
inaccurate filesizes:

http://stackoverflow.com/questions/1392413/calculating-a-directory-size-using-python
'''
# TODO: Manage depth of dirs
# TODO: Filter out by .bcl files and/or include other ones

import os
import sys
import json
import subprocess

def get_dirsize(path="."):
    ''' Gets directory size.
        TODO: Be replaced with a more pythonic way which reports the size correctly.
    '''
    path = path.strip().split('\t')
    out = subprocess.check_output(["du", "-s", path[0]])
    return out.split('\t')[0]

if __name__ == "__main__":
    root = sys.argv[1]
    dirsizes = {}

    for d in os.listdir(root):
    	path = os.path.join(root, d)
        dirsizes[path] = get_dirsize(path)

print json.dumps(dirsizes, sort_keys=True)
