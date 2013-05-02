"""
Halo extension
"""
import os
from cement.core import controller, handler
from scilifelab.pm.core.controller import AbstractBaseController
from scilifelab.utils.halo import run_halo
from scilifelab.utils.misc import query_yes_no

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
        group.add_argument('--setup', help="setup and initialize configuration files", default=False, action="store_true")
        group.add_argument('--config', help="run a given config file", default=None, action="store")
        super(HaloController, self)._setup(app)

    @controller.expose(help="Run halo analysis")
    def run_halo(self):
        if self.app.pargs.setup:
            if not self._check_pargs(["project", "baits", "targets", "target_region"]):
                return
        else:
            if not self._check_pargs(["project"]):
                return
        basedir = os.path.abspath(os.path.join(self.app.controller._meta.root_path, self.app.controller._meta.path_id))
        self.app.log.info("Going to look for samples in {}".format(basedir))
        param_list = run_halo(path=basedir, **vars(self.pargs))
        if self.app.pargs.setup:
            self.app.log.info("Setup configuration files. Rerun command without '--setup' option to run analysis")
            return
        if not len(param_list) > 0:
            self.log.info("No samples found in {}; perhaps you need to add the '--data' option to look in the {} directory".format(self.app.pargs.project, os.path.join(self.app.pargs.project, "data")))
        if len(param_list) > 0 and not query_yes_no("Going to start {} jobs... Are you sure you want to continue?".format(len(param_list)), force=self.pargs.force):
            return
        for param in param_list:
            self.app.cmd.command(param['cl'], **param)
            



