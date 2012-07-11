"""
Pm Archive module

Define a main abstract base controller with arguments common to all
subcommands in the archive group. Add controllers for commands that
require extra arguments.
"""
from cement.core import controller
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

    @controller.expose(hide=True)
    def default(self):
        pass

    @controller.expose(help="List contents of archive folder")
    def ls(self):
        """List contents of archive folder"""
        print "ls"


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


