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

    def sbatch(self,  cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        """sbatch: write cmd_args to sbatch template and submit"""
        print "FIX ME: sbatch: write cmd_args to sbatch template and submit"

    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        ## Is there no easier way to get at --drmaa and --sbatch?!?
        if '--drmaa' in self.app._meta.argv:
            self.drmaa(cmd_args, capture, ignore_error, cwd, **kw)
        elif '--sbatch' in self.app._meta.argv:
            self.sbatch(cmd_args, capture, ignore_error, cwd, **kw)
        else:
            pass

    def drmaa(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        if not self.app.pargs.job_account:
            self.app.log.warn("no job account provided; cannot proceed with drmaa command")
            return
        command = " ".join(cmd_args)
        def runpipe():
            if not os.getenv("DRMAA_LIBRARY_PATH"):
                self.app.log.info("No environment variable DRMAA_LIBRARY_PATH: will not attempt to submit job via DRMAA")
                return
            else:
                import drmaa
            s = drmaa.Session()
            s.initialize()

            jt = s.createJobTemplate()
            jt.remoteCommand = cmd_args[0]
            jt.args = cmd_args[1:]
            jt.jobName = self.app.pargs.jobname
            jt.nativeSpecification = "-A {} -p {} -t {}".format(self.app.pargs.job_account, self.app.pargs.partition, self.app.pargs.time)

            jobid = s.runJob(jt)
            self.app.log.info('Your job has been submitted with id ' + jobid)
    
            s.deleteJobTemplate(jt)
            s.exit()
        return self.dry(command, runpipe)

def add_drmaa_option(app):
    """
    Adds the '--drmaa' argument to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--drmaa', dest='cmd_handler', 
                          action='store_const', help='toggle drmaa command handler', const='drmaa')

def add_sbatch_option(app):
    """
    Adds the '--sbatch' argument to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--sbatch', dest='cmd_handler', 
                          action='store_const', help='toggle sbatch command handler', const = 'sbatch')

def add_shared_distributed_options(app):
    """
    Adds shared distributed arguments to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('-A', '--job_account', type=str,
                          action='store', help='job account', default=None)
    app.args.add_argument('--jobname', type=str,
                          action='store', help='job name', default="pm_distributed")
    app.args.add_argument('-t', '--time',
                          action='store', help='time limit', default="00:10:00")
    app.args.add_argument('--partition', type=str,
                          action='store', help='partition', default="core")

def set_distributed_handler(app):
    """
    Overrides the configured command handler if ``--drmaa`` is passed at the
    command line.
    
    :param app: The application object.
    
    """
    if '--drmaa' in app._meta.argv:
        app._meta.cmd_handler = 'distributed'
        app._setup_cmd_handler()
    elif '--sbatch' in app._meta.argv:
        app._meta.cmd_handler = 'distributed'
        app._setup_cmd_handler()

def load():
    """Called by the framework when the extension is 'loaded'."""
    hook.register('post_setup', add_drmaa_option)
    hook.register('post_setup', add_sbatch_option)
    hook.register('post_setup', add_shared_distributed_options)
    hook.register('pre_run', set_distributed_handler)
    handler.register(DistributedCommandHandler)
