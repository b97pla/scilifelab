"""Pm Archive Module"""

import os
import yaml

from cement.core import controller
from scilifelab.pm.core.controller import AbstractExtendedBaseController
from scilifelab.bcbio.flowcell import *
from scilifelab.lib.archive import flowcell_remove_status

## Main archive controller
class ArchiveController(AbstractExtendedBaseController):
    """
    Functionality for archive management.

    This is the base controller for archive management.
    """
    class Meta:
        """Controller meta-data settings"""
        label = 'archive'
        description = 'Manage archive'

    def _setup(self, base_app):
        super(ArchiveController, self)._setup(base_app)
        base_app.args.add_argument('--as_yaml', action="store_true", default=False, help="list runinfo as yaml file")
        base_app.args.add_argument('-P', '--list-projects', action="store_true", default=False, help="list projects of flowcell")

    def _process_args(self):
        # Set root path for parent class
        self._meta.root_path = self.app.config.get("archive", "root")
        assert os.path.exists(self._meta.root_path), "No such directory {}; check your archive config".format(self._meta.root_path)
        ## Set path_id for parent class
        if self.pargs.flowcell:
            self._meta.path_id = self.pargs.flowcell
        super(ArchiveController, self)._process_args()

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="List runinfo contents")
    def runinfo(self):
        """List runinfo for a given flowcell"""
        if not self._check_pargs(["flowcell"]):
            return
        assert self.config.get("archive", "root"), "archive directory not defined"
        fc = Flowcell(os.path.join(self.config.get("archive", "root"), self.pargs.flowcell, "run_info.yaml"))
        self.log.info("Opening file {}".format(fc.filename))
        if self.pargs.list_projects:
            self.app._output_data['stdout'].write("available projects for flowcell {}:\n\t".format(self.pargs.flowcell) + "\n\t".join(fc.projects()))
            return
        if self.pargs.as_yaml:
            self.app._output_data['stdout'].write(str(fc.as_yaml()))
        else:            
            self.app._output_data['stdout'].write(str(fc))

    @controller.expose(help="Verify flowcells that can be deleted from archive")
    def rm_status(self):
        """This function looks for flowcells that could be deleted
        from archive and returns a list of flowcells with a KEEP/RM
        flag."""
        out_data = flowcell_remove_status(self.app.config.get("archive", "root"), self.app.config.get("production", "swestore"))
        self.app._output_data['stdout'].write(out_data['stdout'].getvalue())
        self.app._output_data['stderr'].write(out_data['stderr'].getvalue())

