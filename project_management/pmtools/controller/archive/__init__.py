"""
Pm Archive module
"""
from cement.core import controller
from pmtools import AbstractBaseController, PmController

## Main archive controller
class ArchiveController(AbstractBaseController):
    """
    Functionality for management of archive folder.
    """
    class Meta:
        label = 'archive'
        stacked_on = 'pm'
        description = 'Manage archive folder'
        interface = controller.IController
        arguments = [
            (['-p', '--project'], dict(help="project"))
            ]

    @controller.expose(hide=True)
    def default(self):
        __doc__


    # @controller.expose(help="Compress files of an archive")
    # def compress(self):
    #     print "compress"


class LsController(controller.CementBaseController):
    class Meta:
        label = 'list'
        stacked_on= 'archive'
        description='List contents of archive folder'
        interface = controller.IController
        arguments = [
            (['-a', '--all'], dict(help="list all")),
            ]

    @controller.expose(hide=True)
    def default(self):
        pass

    # @controller.expose(help="List contents of archive folder")
    # def ls(self):
    #     """List contents of archive folder"""
    #     print "ls"


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


