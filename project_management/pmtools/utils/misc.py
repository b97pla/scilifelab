"""misc"""
import sys
import os
import re

## yes or no: http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input
def query_yes_no(question, default="yes", force=False):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user). The force option simply
    sets the answer to default.
    
    The "answer" return value is one of "yes" or "no".
    
    :param: question
    :param: default
    :param: force
    :returns: yes or no
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        if not force:
            choice = raw_input().lower()
        else:
            choice = "yes"
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                                 "(or 'y' or 'n').\n")

def walk(rootdir):
    """Perform a directory walk"""
    flist = []
    for root, dirs, files in os.walk(rootdir):
        flist = flist + [os.path.join(root, x) for x in files]
    return flist
            
def filtered_walk(rootdir, filter_fn):
    """Perform a filtered directory walk.

    :param rootdir: - root directory
    :param filter_fn: - filtering function
    """
    flist = []
    for root, dirs, files in os.walk(rootdir):
        flist = flist + [os.path.join(root, x) for x in filter(filter_fn, files)]
    return flist
