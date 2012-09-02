"""command"""
import sys

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
        if self.app.pargs.dry_run:
            ##print >> sys.stderr, "(DRY_RUN): " + message + "\n"
            self.app._output_data["stderr"].write("(DRY_RUN): " + message + "\n")
            return
        self.app.log.info(message)
        return func(*args, **kw)
