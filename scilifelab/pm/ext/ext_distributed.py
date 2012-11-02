"""Distributed Extension"""
import re
import os
import sys
import drmaa
import itertools
import argparse

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

    def _save_job_id(self, idfile="JOBID", **job_args):
        """Save jobid to file in working directory"""
        JOBIDFILE = os.path.join(job_args['workingDirectory'], idfile)
        with open(JOBIDFILE, "w") as fh:
            self.app.log.info("Saving jobid {} to file {}".format(self._meta.jobid, JOBIDFILE))
            fh.write(self._meta.jobid)

    def monitor(self, work_dir, idfile="JOBID"):
        """Check for existing job"""
        return self._monitor_job(idfile, **{'workingDirectory':work_dir})

    def _monitor_job(self, idfile="JOBID", **job_args):
        """Check if job is currently being run or in queue. For now,
        the user will manually have to terminate job before proceeding"""
        JOBIDFILE = os.path.join(job_args['workingDirectory'], idfile)
        if not os.path.exists(JOBIDFILE):
            return 
        self.app.log.debug("Will read {} for jobid".format(JOBIDFILE))
        with open(JOBIDFILE) as fh:
            jobid = fh.read()
        ## http://code.google.com/p/drmaa-python/wiki/Tutorial
        decodestatus = {
            drmaa.JobState.UNDETERMINED: 'process status cannot be determined',
            drmaa.JobState.QUEUED_ACTIVE: 'job is queued and active',
            drmaa.JobState.SYSTEM_ON_HOLD: 'job is queued and in system hold',
            drmaa.JobState.USER_ON_HOLD: 'job is queued and in user hold',
            drmaa.JobState.USER_SYSTEM_ON_HOLD: 'job is queued and in user and system hold',
            drmaa.JobState.RUNNING: 'job is running',
            drmaa.JobState.SYSTEM_SUSPENDED: 'job is system suspended',
            drmaa.JobState.USER_SUSPENDED: 'job is user suspended',
            drmaa.JobState.DONE: 'job finished normally',
            drmaa.JobState.FAILED: 'job finished, but failed',
            }
        s = drmaa.Session()
        s.initialize()
        try:
            status = s.jobStatus(str(jobid))
            self.app.log.debug("Getting status for jobid {}".format(jobid))
            self.app.log.info("{}".format(decodestatus[status]))
            if status in [drmaa.JobState.QUEUED_ACTIVE, drmaa.JobState.RUNNING, drmaa.JobState.UNDETERMINED]:
                self.app.log.warn("{}; please terminate job before proceeding".format(decodestatus[status]))
                return True
        except drmaa.errors.InternalException:
            self.app.log.warn("No such jobid {}".format(jobid))
            pass
        s.exit()
        return


    def drmaa(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        if self.app.pargs.partition == "node" and self.app.pargs.max_node_jobs < self._meta.n_submitted_jobs:
            self.app.log.info("number of submitted jobs larger than maximum number of allowed node jobs; not submitting job")
            return
        self._meta.n_submitted_jobs = self._meta.n_submitted_jobs + 1
        if not self._check_args(**kw):
            self.app.log.warn("missing argument; cannot proceed with drmaa command. Make sure you provide time, account, partition, and jobname")
            return
        if kw.get('platform_args', None):
            platform_args = opt_to_dict(kw['platform_args'])
        else:
            platform_args = opt_to_dict([])
        kw.update(**vars(self.app.pargs))
        job_args = make_job_template_args(platform_args, **kw)
        if kw.get('monitorJob', False):
            if self._monitor_job(**job_args):
                self.app.log.info("exiting from {}".format(__name__) )
                return

        command = " ".join(cmd_args)
        def runpipe():
            s = drmaa.Session()
            s.initialize()
            jt = s.createJobTemplate()
            jt.remoteCommand = cmd_args[0]
            jt.args = cmd_args[1:]
            jt.jobName = job_args['jobname']
            if os.path.isdir(job_args['outputPath']):
                jt.outputPath = ":" + drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.join(os.path.relpath(job_args['outputPath'], os.getenv("HOME")), jt.jobName + "-drmaa.log")
            else:
             jt.outputPath = ":" + drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.join(os.path.relpath(job_args['outputPath'], os.getenv("HOME")))
            jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.relpath(job_args['workingDirectory'], os.getenv("HOME"))

            jt.nativeSpecification = "-t {time} -p {partition} -A {account} {extra}".format(**job_args)
            if kw.get('email', None):
                jt.email=[kw.get('email')]
            self.app.log.info("Submitting job with native specification {}".format(jt.nativeSpecification))
            self._meta.jobid = s.runJob(jt)
            self.app.log.info('Your job has been submitted with id ' + self._meta.jobid)
            if kw.get('saveJobId', False):
                self._save_job_id(**job_args)
            s.deleteJobTemplate(jt)
            s.exit()
            
        return self.dry(command, runpipe)

def make_job_template_args(opt_d, **kw):
    """Given a dictionary of arguments, update with kw dict that holds arguments passed to argv.

    :param opt_d: dictionary of option key/value pairs
    :param kw: dictionary of program arguments

    :returns: dictionary of job arguments
    """
    job_args = {}
    job_args['jobname'] = kw.get('jobname', None) or opt_d.get('-J', None) or  opt_d.get('--job-name', None)
    job_args['time'] = kw.get('time', None) or opt_d.get('-t', None) or  opt_d.get('--time', None)
    job_args['partition'] = kw.get('partition', None) or opt_d.get('-p', None) or  opt_d.get('--partition', None)
    job_args['account'] = kw.get('account', None) or opt_d.get('-A', None) or  opt_d.get('--account', None)
    job_args['outputPath'] = kw.get('outputPath', None) or opt_d.get('-o', os.curdir)
    job_args['workingDirectory'] = kw.get('workingDirectory', None) or opt_d.get('-D', None) 
    job_args['email'] = kw.get('email', None) or opt_d.get('--mail-user', None) 
    invalid_keys = ["--mail-user", "--mail-type", "-o", "--output", "-D", "--workdir", "-J", "--job-name", "-p", "--partition", "-t", "--time", "-A", "--account"]
    extra_keys = [x for x in opt_d.keys() if x not in invalid_keys]
    extra_args = ["{}={}".format(x, opt_d[x]) if x.startswith("--") else "{} {}".format(x, opt_d[x]) for x in extra_keys]
    job_args['extra'] = kw.get('extra_args', None) or extra_args
    job_args['extra'] = " ".join(job_args['extra'])
    return job_args

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
    group.add_argument('--extra_args', type=str,  nargs=argparse.REMAINDER,
                          action='store', help='extra arguments to pass to drmaa native specification. NOTE: must be supplied last since it uses remaining part of argument list', default=None)
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
