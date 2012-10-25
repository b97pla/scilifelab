"""Distributed Extension"""
import re
import os
import sys
import drmaa
import itertools

from cement.core import backend, handler, hook

from scilifelab.pm.core import command

LOG = backend.minimal_logger(__name__)

class DistributedCommandHandler(command.CommandHandler):
    """ 
    This class is an implementation of the :ref:`ICommand
    <scilifelab.pm.core.command>` interface.
    """    

    class Meta:
        """Handler meta-data"""
        
        interface = command.ICommand
        """The interface that this class implements."""

        label = 'distributed'
        """The string identifier of this handler."""

        n_submitted_jobs = 0
        """The number of submitted jobs"""

        jobid = None
        """The submitted jobid"""

        platform_args = None
        """Platform specific arguments"""

    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        ## Is there no easier way to get at --drmaa?!?
        if '--drmaa' in self.app._meta.argv:
            self.drmaa(cmd_args, capture, ignore_error, cwd, **kw)
        else:
            pass

    def _check_args(self, **kw):
        pargs = kw.get('platform_args', [])
        if not self.app.pargs.account and "-A" not in pargs and "--account" not in pargs:
            return False
        if not self.app.pargs.jobname and "-J" not in pargs and "--jobname" not in pargs:
            return False
        if not self.app.pargs.partition and "-p" not in pargs and "--partition" not in pargs:
            return False
        if not self.app.pargs.time and "-t" not in pargs and "--time" not in pargs:
            return False
        return True

    def drmaa(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        if self.app.pargs.partition == "node" and self.app.pargs.max_node_jobs < self._meta.n_submitted_jobs:
            self.app.log.info("number of submitted jobs larger than maximum number of allowed node jobs; not submitting job")
            return
        self._meta.n_submitted_jobs = self._meta.n_submitted_jobs + 1
        if not self._check_args(**kw):
            self.app.log.warn("missing argument; cannot proceed with drmaa command. Make sure you provide time, account, partition, and jobname")
            return
        command = " ".join(cmd_args)
        def runpipe():
            s = drmaa.Session()
            s.initialize()
            jt = s.createJobTemplate()
            jt.remoteCommand = cmd_args[0]
            jt.args = cmd_args[1:]
            if kw['platform_args']:
                platform_args = opt_to_dict(kw['platform_args'])
            else:
                platform_args = opt_to_dict([])
            opt_d = make_job_template_args(platform_args, **vars(self.app.pargs))
            jt.outputPath = ":" + drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.relpath(opt_d['outputPath'], os.getenv("HOME"))
            jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.relpath(opt_d['workingDirectory'], os.getenv("HOME"))
            jt.jobName = opt_d['jobname']
            jt.nativeSpecification = "-t {time} -p {partition} -A {account}".format(**opt_d)
            self.app.log.info("Submitting job with native specification {}".format(jt.nativeSpecification))
            self._meta.jobid = s.runJob(jt)
            self.app.log.info('Your job has been submitted with id ' + self._meta.jobid)
            s.deleteJobTemplate(jt)
            s.exit()
            
        return self.dry(command, runpipe)


def opt_to_dict(opts):
    """Transform option list to a dictionary.

    :param opts: option list
    
    :returns: option dictionary
    """
    if isinstance(opts, dict):
        return
    args = list(itertools.chain.from_iterable([x.split("=") for x in opts]))
    opt_d = {k: True if v.startswith('-') else v
             for k,v in zip(args, args[1:]+["--"]) if k.startswith('-')}
    return opt_d

def make_job_template_args(opt_d, **kw):
    """Given a dictionary of arguments, update with kw dict that holds arguments passed to argv.

    :param opt_d: dictionary of option key/value pairs
    :param kw: dictionary of keywords
    """
    kw['jobname'] = kw.get('jobname', None) or opt_d.get('-J', None) or  opt_d.get('--job-name', None)
    kw['time'] = kw.get('time', None) or opt_d.get('-t', None) or  opt_d.get('--time', None)
    kw['partition'] = kw.get('partition', None) or opt_d.get('-p', None) or  opt_d.get('--partition', None)
    kw['account'] = kw.get('account', None) or opt_d.get('-A', None) or  opt_d.get('--account', None)
    kw['outputPath'] = kw.get('outputPath', os.curdir) or opt_d.get('-o', None)
    kw['workingDirectory'] = kw.get('workingDirectory', None) or opt_d.get('-D', None) 
    return kw

def add_drmaa_option(app):
    """
    Adds the '--drmaa' argument to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--drmaa', dest='cmd_handler', 
                          action='store_const', help='toggle drmaa command handler', const='drmaa')

def add_shared_distributed_options(app):
    """
    Adds shared distributed arguments to the argument object.
    
    :param app: The application object.
    
    """
    group = app.args.add_argument_group('distributed', 'Options for distributed execution.')
    group.add_argument('-A', '--account', type=str,
                          action='store', help='job account', default=None)
    group.add_argument('--jobname', type=str,
                          action='store', help='job name', default=None)
    group.add_argument('-t', '--time',
                          action='store', help='time limit', default=None)
    group.add_argument('--partition', type=str,
                          action='store', help='partition (node, core or devel)', default=None)
    group.add_argument('--max_node_jobs', type=int, default=10,
                          action='store', help='maximum number of node jobs (default 10)')

def set_distributed_handler(app):
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
    if not os.getenv("DRMAA_LIBRARY_PATH"):
        self.app.log.warn("No environment variable $DRMAA_LIBRARY_PATH: loading {} failed".format(__name__))
        return
    hook.register('post_setup', add_drmaa_option)
    hook.register('post_setup', add_shared_distributed_options)
    hook.register('pre_run', set_distributed_handler)
    handler.register(DistributedCommandHandler)
