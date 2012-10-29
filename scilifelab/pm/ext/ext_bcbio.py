"""bcbio extension"""
import os
import re

from cement.core import backend, controller, handler, hook
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.misc import query_yes_no, filtered_walk

class BcbioBaseController(AbstractBaseController):
    """
    Basic functionality for dealing with bcbio based analyses.

    """

    class Meta:
        label = 'bcbio_base'
        description = 'Extension for dealing with bcbio'

    def _setup(self, app):
        group = app.args.add_argument_group('Bcbio argument group', 'Options for bcbio')
        group.add_argument('-M', '--merge_analysis', action='store_true',
                       help='Run merge analysis. If a sample has data from more than one run, a \'total\' directory will be setup.', default=False)
        
        super(BcbioBaseController, self)._setup(app)

    @controller.expose(help="Run bcbb pipeline")
    def run2(self):
        if not self._check_pargs(["project"]):
            return
        print self._help_text


class BcbioProjectController(BcbioBaseController):
    """
    Basic functionality for dealing with bcbio based analyses in projects.
    """

    class Meta:
        label = 'bcbio_project'
        stacked_on = 'project'
        

class BcbioProductionController(BcbioBaseController):
    """
    Basic functionality for dealing with bcbio based analyses in production.
    """

    class Meta:
        label = 'bcbio_production'
        stacked_on = 'production'


def load():
    """Called by the framework when the extension is 'loaded'."""
    handler.register(BcbioProjectController)
    handler.register(BcbioProductionController)
