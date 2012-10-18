"""Scilifelab dry module"""
import os

import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)
from cement.core import backend
LOG = backend.minimal_logger(__name__)

def dry(message, func, dry_run=True, *args, **kw):
    """Wrapper that runs a function (runpipe) if flag dry_run isn't set, otherwise returns function call as string
    
    :param message: message describing function call
    :param func: function to call
    :param *args: positional arguments to pass to function
    :param **kw: keyword arguments to pass to function
    """
    if dry_run:
        LOG.debug("(DRY_RUN): " + message + "\n")
        return "(DRY_RUN): " + message + "\n"
    return func(*args, **kw)

def dry_unlink(fh, dry_run=True):
    """Wrapper for unlinking a file.
    
    :param fh: file name
    :param dry_run: flag to set dryness
    """
    def runpipe():
        if fh is None:
            return
        if not os.path.exists(fh):
            LOG.warn("not going to remove non-existant file {}".format(fh))
            return
        os.unlink(fh)
    return dry("removing file {}".format(fh), runpipe, dry_run)

def dry_write(fn, data=None, dry_run=True):
    """Wrapper for writing data to a file.
    
    :param fn: file name <str>
    :param data: data structure to write as <str>
    :param dry_run: flag to set dryness
    """
    def runpipe():
        if fn is None:
            return
        if os.path.exists(fn):
            LOG.warn("not overwriting existing file {}".format(fn))
            return
        with open (fn, "w") as fh:
            fh.write(data)
    return dry("writing data to file {}".format(fn), runpipe, dry_run)
