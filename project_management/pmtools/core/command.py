"""
pmtools command core module.

"""
import os
import sys
import shutil

from cement.core import interface, handler

def cmd_interface_validator(cls, obj):
    members = [
        '_setup',
        'command',
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
            ##print >> sys.stderr, "(DRY_RUN): " + message + "\n"
            self.app._output_data["stderr"].write("(DRY_RUN): " + message + "\n")
            return
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

    def transfer_file(self, src, tgt):
        """Wrapper for transferring files with move or copy operation.

        :param src: source destination
        :param tgt: target destination
        """
        if self.app.pargs.move:
            deliver_fn = shutil.move
        else:
            deliver_fn = shutil.copyfile
        def runpipe():
            if src is None:
                return
            if os.path.exists(tgt):
                self.app.log.warn("{} already exists: not doing anything!".format(tgt))
                return
            deliver_fn(src, tgt)
        return self.dry("{} file {} to {}".format(deliver_fn.__name__, src, tgt), runpipe) 

    def write(self, fn, data):
        """Wrapper for writing data to a file.

        :param fn: file name <str>
        :param data: data structure to write as <str>
        """
        def runpipe():
            if data is None:
                return
            with open (fn, "w") as fh:
                fh.write(data)
        return self.dry("writing data to file {}".format(fn), runpipe)
