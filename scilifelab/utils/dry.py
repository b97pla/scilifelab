"""Scilifelab dry module"""
import os
import shutil

import scilifelab.log
LOG = scilifelab.log.minimal_logger(__name__)

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

def dry_unlink(fn, dry_run=True):
    """Wrapper for unlinking a file.
    
    :param fn: file name
    :param dry_run: flag to set dryness
    """
    def runpipe():
        if fn is None:
            return
        if not os.path.exists(fn):
            LOG.warn("not going to remove non-existant file {}".format(fn))
            return
        os.unlink(fn)
    return dry("removing file {}".format(fn), runpipe, dry_run)

def dry_rmdir(dname, dry_run=True):
    """Wrapper for removing a dir.
    
    :param dname: dir name
    :param dry_run: flag to set dryness
    """
    def runpipe():
        if dname is None:
            return
        if not os.path.exists(dname):
            LOG.warn("not going to remove non-existant directory {}".format(dname))
            return
        if os.path.isdir(dname):
            try:
                os.rmdir(dname)
            except:
                pass
    return dry("removing directory {}".format(dname), runpipe, dry_run)

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

def dry_backup(fn, ext=".bak", dry_run=True):
    """Wrapper for making a backup copy of a file.

    :param fn: file name <str>
    :param ext: backup extension
    :param dry_run: dry run flag
    """
    fn_bak = "{}{}".format(fn, ext)
    def runpipe():
        if fn is None:
            return
        if os.path.exists(fn_bak):
            LOG.warn("not overwriting existing file {}".format(fn_bak))
            return 
        shutil.copyfile(fn, fn_bak)
    return dry("making backup of file {} as ".format(fn, fn_bak), runpipe, dry_run)
    
