"""pmtools command core module."""
import os
import sys
import shutil

from cement.core import interface, handler

def cmd_interface_validator(cls, obj):
    members = [
        '_setup',
        'command',
        'monitor',
        ]
    interface.validate(ICommand, obj, members)

class ICommand(interface.Interface):

    class IMeta:
        """Interface meta-data"""

        label = 'command'
        """The string identifier of the interface."""

        validator = cmd_interface_validator
        """The interface validator function."""

    Meta = interface.Attribute('Handler Meta-data')

    def _setup(app_obj):
        """
        The _setup function is called during application initialization and
        must 'setup' the handler object making it ready for the framework
        or the application to make further calls to it.
        
        :param app_obj: The application object. 
                                
        """

    def monitor(work_dir, idfile=None):
        """
        Check for process/job id file.
        
        :param work_dir: working directory
        :param idfile: process/job id file
        """

    def command(cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        """
        Run the command interface function

        Required arguments:

          cmd_args
            Command arguments <list>
          capture
            Capture output
          ignore_error
            Ignore errors
          cwd
            Working directory
          

        Optional arguments:

          kw
            keyword arguments
           
        """
        def runpipe():
            """
            The actual function to run
            """

class CommandHandler(handler.CementBaseHandler):
    """
    Base class that all Command Handlers should sub-class from.

    """
    class Meta:
        """
        Handler meta-data
        """

        label = None
        """The string identifier of this handler."""

        interface = ICommand
        """The interface that this class implements."""

    def __init__(self, *args, **kw):
        super(CommandHandler, self).__init__(*args, **kw)

    ## Taken from paver.easy
    ## FIXME: add time stamp (better: make DRY_RUN a log level that only prints to console, for instance by using the interface ILog)
    def dry(self, message, func, *args, **kw):
        """Wrapper that runs a function (runpipe) if flag dry_run isn't set, otherwise returns function call as string

        :param message: message describing function call
        :param func: function to call
        :param *args: positional arguments to pass to function
        :param *kw: keyword arguments to pass to function
        """
        if self.app.pargs.dry_run:
            self.app._output_data["stderr"].write("(DRY_RUN): " + message + "\n")
            return
        if self.app.pargs.verbose:
            self.app.log.info(message)
        return func(*args, **kw)

    def safe_makedir(self, dname):
        """Make a directory if it doesn't exist.
        
        :param dname: directory to create
        """
        def runpipe():
            if not os.path.exists(dname):
                try:
                    os.makedirs(dname)
                except OSError:
                    if not os.path.isdir(dname):
                        raise
            else:
                self.app.log.info("Directory %s already exists" % dname)
            return dname
        return self.dry("Make directory %s" % dname, runpipe)

    def _rsync(self, src, tgt):
        """Wrapper for running rsync.

        :param src: source destination
        :param tgt: target destination
        """
        opts = "-av"
        if os.path.isdir(src):
            src = os.path.join(src) + os.sep
        if os.path.isdir(tgt):
            tgt = os.path.join(tgt) + os.sep
        cl = ["rsync {} {} {}".format(opts, src, tgt)]
        out = self.app.cmd.command(cl, **{'shell':True})
        

    def transfer_file(self, src, tgt):
        """Wrapper for transferring files with move or copy operation.

        :param src: source destination
        :param tgt: target destination
        """
        if self.app.pargs.move:
            deliver_fn = shutil.move
        elif self.app.pargs.copy:
            deliver_fn = shutil.copyfile
        elif self.app.pargs.rsync:
            deliver_fn = self._rsync
        else:
            pass
        def runpipe():
            if src is None:
                return
            if not self.app.pargs.rsync and os.path.exists(tgt):
                self.app.log.warn("{} already exists: not doing anything!".format(tgt))
                return
            deliver_fn(src, tgt)
        return self.dry("{} file {} to {}".format(deliver_fn.__name__, src, tgt), runpipe) 

    def write(self, fn, data=None, overwrite=False):
        """Wrapper for writing data to a file.

        :param fn: file name <str>
        :param data: data structure to write as <str>
        """
        def runpipe():
            if fn is None:
                return
            if os.path.exists(fn) and not overwrite:
                self.app.log.warn("not overwriting existing file {}".format(fn))
                return
            with open (fn, "w") as fh:
                fh.write(data)
        return self.dry("writing data to file {}".format(fn), runpipe)
    
    def safe_touchfile(self, fname):
        """Touch a non-existing file
        
        :param fname: file name
        """
        return self.write(fname,data="",overwrite=False)
                
    def safe_unlink(self, fh):
        """Wrapper for unlinking a file.

        :param fh: file name
        """
        def runpipe():
            if fh is None:
                return
            if not os.path.exists(fh):
                self.app.log.warn("not going to remove non-existant file {}".format(fh))
                return
            os.unlink(fh)
        return self.dry("removing file {}".format(fh), runpipe)

    def safe_rmdir(self, d):
        """Wrapper for removing directories. Directory structure
        should be empty of files, possibly containing subdirectories.

        :param d: directory
        """
        def runpipe():
            if d is None:
                return
            if not os.path.exists(d):
                self.app.log.warn("not going to remove non-existant directory {}".format(d))
                return
            dirlist = []
            for root, dirs, files in os.walk(d):
                dirlist.extend([os.path.join(root, x) for x in dirs])
            for x in dirlist:
                try:
                    os.removedirs(x)
                except:
                    pass
        return self.dry("removing directory {}".format(d), runpipe)

    def rmtree(self, d):
        """Wrapper for recursively removing a directory structure, possibly containing
        files and subdirectories
        
        :param d: directory
        """
        def runpipe():
            if d is None:
                return
            if not os.path.exists(d):
                self.app.log.warn("not going to remove non-existant directory {}".format(d))
                return
            self.app.log.debug("removing directory {} and everything beneath it".format(d))
            try:
                shutil.rmtree(d)
            except e:
                self.app.log.warn(e)
        return self.dry("removing directory {} and everything beneath it".format(d),runpipe)
    
    def link(self, src, tgt):
        """Wrapper for making links.

        :param src: source link
        :param tgt: target link
        """
        def runpipe():
            if not os.path.exists(tgt):
                try:
                    os.symlink(src, tgt)
                except:
                    self.app.log.warn("Couldn't create link {} -> {}".format(tgt, src))
                    pass
        return self.dry("creating link {} -> {}".format(tgt, src), runpipe)

    def chown(self, fname, uid=-1, gid=-1):
        """Wrapper for changing ownership of a file
        
        :param fname: file name
        :param uid: uid to set (-1 for unchanged)
        :param gid: gid to set (-1 for unchanged)
        """
        def runpipe():
            if not os.path.exists(fname):
                self.app.log.warn("not changing ownership of non-existant file {}".format(fname))
                return
            os.chown(fname,uid,gid)
        return self.dry("changing ownership of file {}".format(fname), runpipe)

    def chmod(self, fname, mode):
        """Wrapper for changing the mode of a file
        
        :param fname: file name
        :param mode: numeric mode to set
        """
        def runpipe():
            if not os.path.exists(fname):
                self.app.log.warn("not changing mode of non-existant file {}".format(fname))
                return
            os.chmod(fname,mode)
        return self.dry("changing mode of file {}".format(fname), runpipe)

    def md5sum(self, fname):
        """Calculate the md5sum of a file and write the output to a file or, if supplied, an output pipe
        """
        import scilifelab.utils.misc
        def runpipe():
            if not os.path.exists(fname):
                self.app.log.warn("not calculating md5sum of non-existant file {}".format(fname))
                return
            self.app.log.debug("Calculating md5sum of file {}".format(fname))
            md5 = scilifelab.utils.misc.md5sum(fname)
            md5file = "{}.md5".format(fname)
            self.app.log.debug("Writing md5sum to file {}".format(md5file))
            self.write(md5file,"{}  {}".format(md5,os.path.basename(fname)),True)
            return md5file
        return self.dry("calculating md5sum of {}".format(fname), runpipe)
    
    def verify_md5sum(self, md5file):
        """Verify the md5sums and files given in the supplied md5file
        """
        import scilifelab.utils.misc
        def runpipe():
            if not os.path.exists(md5file):
                self.app.log.warn("not verifying md5sums in non-existant file {}".format(md5file))
                return False
            passed = True
            self.app.log.debug("Verifying md5sums in file {}".format(md5file))
            with open(md5file) as fh:
                for line in fh:
                    line = line.strip()
                    if len(line) == 0:
                        continue
                    self.app.log.debug("Processing line: {}".format(line))
                    pcs = line.split()
                    if not len(pcs) == 2:
                        self.app.log.warn("malformed line: {} in {}".format(line,md5file))
                        continue
                    fpath = os.path.join(os.path.dirname(md5file),pcs[1])
                    self.app.log.debug("Calculating md5sum of file {}. Expecting {}".format(fpath,pcs[0]))
                    md5 = scilifelab.utils.misc.md5sum(fpath)
                    self.app.log.debug("Calculated md5sum is {}".format(md5))
                    if md5 == pcs[0]:
                        self.app.log.info("{}: OK".format(pcs[1]))
                    else: 
                        self.app.log.warn("{}: FAILED".format(pcs[1]))
                        passed = False
            return passed
        return self.dry("verifying md5sums in {}".format(md5file), runpipe)
    