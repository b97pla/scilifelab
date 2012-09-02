"""Distributed Extension"""

import os
import sys

from cement.core import backend, handler, hook

from pmtools.core import command

LOG = backend.minimal_logger(__name__)

class DistributedCommandHandler(command.CommandHandler):
    """ 
    This class is an implementation of the :ref:`ICommand
    <pmtools.core.command>` interface.
    """    

    class Meta:
        """Handler meta-data"""
        
        interface = command.ICommand
        """The interface that this class implements."""

        label = 'distributed'
        """The string identifier of this handler."""


    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        if not os.getenv("DRMAA_LIBRARY_PATH"):
            self.app.log.info("No environment variable DRMAA_LIBRARY_PATH: will not attempt to submit job via DRMAA")
            return
        else:
            import drmaa

        if self.pargs.node:
            partition = "node"
        if not self.app.pargs.uppmax_project:
            self.app.log.warn("no uppmax id provided; cannot proceed with drmaa command")
            sys.exit()
        command = " ".join(cmd_args)
        def runpipe():
            print "in runpipe"
            s = drmaa.Session()
            s.initialize()
            print s
            jt = s.createJobTemplate()
            jt.remoteCommand = cmd_args[0]
            jt.args = cmd_args[1:]

            # # TODO: job name is always (null), must fix slurm_drmaa C library and its
            # # custom parsing (substitute "slurmdrmaa_parse_native"
            # # for GNU GetOpt on slurm_drmaa/util.c)
            jt.job_name = jobname
            jt.nativeSpecification = "-A a2010002 -p devel"# % (self.pargs.uppmax_project, partition)#, str(self.pargs.sbatch_time))

            print jt
            jobid = s.runJob(jt)
            print jobid
            self.log.info('Your job has been submitted with id ' + jobid)
    
            s.deleteJobTemplate(jt)
            s.exit()
            print "Exiting runpipe"
        return self.dry(command, runpipe)
        ##return self._not_implemented("Implement drmaa code as soon as drmaa library fixed!\nSee AbstractBaseController.drmaa function")


def add_drmaa_option(app):
    """
    Adds the '--drmaa' argument to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--drmaa', dest='cmd_handler', 
        action='store_const', help='toggle command handler', const='distributed')

def set_drmaa_handler(app):
    """
    Overrides the configured command handler if ``--drmaa`` is passed at the
    command line.
    
    :param app: The application object.
    
    """
    if '--drmaa' in app._meta.argv:
        app._meta.cmd_handler = 'distributed'
        app._setup_cmd_handler()

def load():
    """Called by the framework when the extension is 'loaded'."""
    hook.register('post_setup', add_drmaa_option)
    hook.register('pre_run', set_drmaa_handler)
    handler.register(DistributedCommandHandler)
