"""
Halo extension
"""
import os
from cement.core import controller, handler
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.halo import run_halo

import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)

class HaloController(AbstractBaseController):
    """
    This class is an implementation of the :ref:`IController
    <scilifelab.pm.core.controller>` interface.

    Functionality for dealing with HaloPlex analyses.
    """
    class Meta:
        label = 'halo'
        description = 'Extension for dealing with halo data'

    def _setup(self, app):
        group = app.args.add_argument_group('HaloPlex argument group', 'Options for halo analyses')
        group.add_argument('--batch_size', help="set the batch size (number of samples run per batch)", default=8, type=int, action="store")
        group.add_argument('--target_region', help="set the target region (must be a bed file)", default=None, type=str, action="store")
        super(HaloController, self)._setup(app)

    @controller.expose(help="Run halo analysis")
    def run_halo(self):
        if not self._check_pargs(["project", "baits", "targets", "target_region"]):
            return
        basedir = os.path.abspath(os.path.join(self.app.controller._meta.root_path, self.app.controller._meta.path_id))
        self.app.log.info("Going to look for samples in {}".format(basedir))
        out_data = run_halo(path=basedir, **vars(self.pargs))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())



