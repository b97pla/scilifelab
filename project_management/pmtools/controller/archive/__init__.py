"""
Pm Archive module

Define a main abstract base controller with arguments common to all
subcommands in the archive group. Add controllers for commands that
require extra arguments.
"""
import re
from cement.core import controller
from cement.utils.shell import *
from pmtools import AbstractBaseController, SubSubController

## Main archive controller
class ArchiveController(AbstractBaseController):
    """
    Functionality for management of archive folder.
    """
    class Meta:
        label = 'archive'
        stacked_on = None
        description = 'Manage archive folder'
        interface = controller.IController
        arguments = [
            (['-f', '--flowcell'], dict(help="Flowcell id")),
            ]

    @controller.expose(hide=True)
    def default(self):
        __doc__


    @controller.expose(help="Compress files of an archive")
    def compress(self):
        print "compress"


class LsController(SubSubController):
    class Meta:
        label = 'list'
        stacked_on= 'archive'
        description='List controller'
        interface = controller.IController
        arguments = [
            (['-a', '--all'], dict(help="list all")),
            ]

    def _setup(self, app_obj):
        # shortcuts
        super(SubSubController, self)._setup(app_obj)
        # Compile ignore regexps
        self.reignore = re.compile(self.config.get("config", "ignore").replace("\n", "|"))

    def _filtered_ls(self, out):
        """Filter output"""
        def ignore(line):
            return self.reignore.match(line) == None
        return filter(ignore, out)
    
    @controller.expose(hide=True)
    def default(self):
        pass

    @controller.expose(help="List contents of archive folder")
    def ls(self):
        """List contents of archive folder"""
        (out, err, code) = exec_cmd(["ls",  self.app.config.get("config", "archive")])
        if code == 0:
            ## FIXME: use output formatter for stuff like this
            print "\n".join(self._filtered_ls(out.splitlines()))
        else:
            self.log.warn(err)


class HelloController(controller.CementBaseController):
    class Meta:
        label = 'helloctrl'
        stacked_on= 'archive'
        description='List contents of archive folder'
        interface = controller.IController

    @controller.expose(hide=True)
    def default(self):
        pass

    # @controller.expose(help="List contents of archive folder")
    # def hello(self):
    #     """List contents of archive folder"""
    #     print "ls"


