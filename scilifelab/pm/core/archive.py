"""Pm Archive Module"""

import os
import yaml
import sys

from cement.core import controller
from scilifelab.pm.core.controller import AbstractExtendedBaseController
from scilifelab.bcbio.flowcell import *
from scilifelab.lib.archive import *

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
        base_app.args.add_argument('--package-run', action="store_true", default=False, help="package a run in preparation for archiving to swestore")
        base_app.args.add_argument('--clean-from-staging', action="store_true", default=False, help="Removes the uncompressed run from staging if archiving was OK")
        base_app.args.add_argument('--check-finished', action="store_true", default=False, help="Whether to check if a run has finished transfer before packing it")
        base_app.args.add_argument('--excludes', action="store", default=None, help="a file containing file and directory name patterns to exclude from tarball")
        base_app.args.add_argument('--workdir', action="store", default=None, help="the folder to create tarballs in (default is parent of the run folder)")
        base_app.args.add_argument('--force-overwrite', action="store_true", default=False, help="will replace any existing package")
        base_app.args.add_argument('--tarball', action="store", default=None, help="tarball package of a run to send to swestore")
        base_app.args.add_argument('--clean', action="store_true", default=False, help="remove the flowcell run folder after successful packaging")
        base_app.args.add_argument('--remote-upload', action="store_true", default=False, help="upload the tarball to the remote host specified in config file")
        base_app.args.add_argument('--remote-user', action="store", default=None, help="remote user to use for remote upload")
        base_app.args.add_argument('--remote-host', action="store", default=None, help="remote host to use for remote upload")
        base_app.args.add_argument('--remote-path', action="store", default=None, help="remote path to use for remote upload")
        base_app.args.add_argument('--send-to-swestore', action="store_true", default=False, help="send a packaged run to swestore using irods")
        base_app.args.add_argument('--clean-swestore', action="store_true", default=False, help="Clean the tarball after successfuly archiving in swestore")
        base_app.args.add_argument('--swestore-path', action="store", default=None, help="the path to the project's folder on the swestore area")
        base_app.args.add_argument('--remote-swestore', action="store_true", default=False, help="run the swestore archiving script on the remote host instead of locally")
        base_app.args.add_argument('--log-to-db', action="store_true", default=False, help="log the swestore archiving progress to db")


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

    @controller.expose(help="Remove a run folder from archive")
    def rm(self):
        """Remove a specified flowcell run folder from archive
        """

        # We require a flowcell argument
        if not self._check_pargs(["flowcell"]):
            return

        rm_run(self,self.config.get('archive','root'), flowcell=self.pargs.flowcell)

    @controller.expose(help="Archive a run in swestore")
    def swestore(self):
        """This function is the entry point for tasks having to do with packaging and sending runs to swestore
        """
        db_info = self.app.config.get_section_dict('db')
        f_conn = FlowcellRunMetricsConnection(username=db_info.get('user'),
                                              password=db_info.get('password'),
                                              url=db_info.get('url'))
        # Create a tarball out of the run folder
        if self.pargs.package_run:

            # We require a flowcell argument
            if not self._check_pargs(["flowcell"]):
                return

            self.pargs.tarball = package_run(self,self.config.get('archive','swestore_staging'), **vars(self.pargs))
            if not self.pargs.tarball:
                self.log.error("No tarball was created, exiting")
                return
            if self.pargs.clean:
                rm_run(self,self.config.get('archive','root'), flowcell=self.pargs.flowcell)

            if self.pargs.clean_from_staging:
                #Check that the run has been archived on the NAS before removing it, otherwise it will keep synching
                if self.pargs.flowcell in f_conn.get_storage_status('NAS_nosync').keys():
                    rm_run(self,self.config.get('archive','swestore_staging'), flowcell=self.pargs.flowcell)
                else:
                    self.log.warn("Run storage status is not NAS_nosync, not removing run from swestore_stage!")

        if not self.pargs.tarball:
            self.log.error("Required argument --tarball was not specified")
            return

        if not os.path.exists(self.pargs.tarball):
            self.log.error("Tarball {} does not exist".format(self.pargs.tarball))
            return

        # Upload a tarball to a remote host
        if self.pargs.remote_upload:
            result = upload_tarball(self,
                                    **dict(self.config.get_section_dict('archive').items() + vars(self.pargs).items()))
            if not result:
                return
            if self.pargs.clean:
                rm_tarball(self,tarball=self.pargs.tarball)

        # Send the tarball to Swestore using irods
        if self.pargs.send_to_swestore:
            result = send_to_swestore(self,**dict(self.config.get_section_dict('archive').items() + vars(self.pargs).items()))
            if not result:
                # If archiving failed, we need to give a non-zero exit code in order for a remote instance to detect the failure
                sys.exit(1)
            if self.pargs.clean:
                rm_tarball(self,tarball=self.pargs.tarball)
            #Set the run as archived in StatusDB
            fc_db_id = f_conn.id_view.get(self.pargs.flowcell)
            f_conn.set_storage_status(fc_db_id, 'swestore_archived')
            # Log to statusdb
            if self.pargs.log_to_db:
                # implement this
                raise NotImplementedError("logging to db functionality not implemented")

