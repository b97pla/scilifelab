"""Distributed Extension"""
import re
import os
import sys
try:
    import drmaa
except:
    pass
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

        batch_command = []
        """Batch command array"""

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
        if not self.app.pargs.jobname and "-J" not in pargs and "--job-name" not in pargs:
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
        if kw.get('platform_args', None):
            platform_args = opt_to_dict(kw['platform_args'])
        else:
            platform_args = opt_to_dict([])
        kw.update(**vars(self.app.pargs))
        job_args = make_job_template_args(platform_args, **kw)
        if not self._check_args(**kw):
            self.app.log.warn("missing argument; cannot proceed with drmaa command. Make sure you provide time, account, partition, and jobname")
            return

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
            if os.path.isdir(job_args['errorPath']):
                jt.errorPath = ":" + drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.join(os.path.relpath(job_args['errorPath'], os.getenv("HOME")), jt.jobName + "-drmaa.err")
            else:
                jt.errorPath = ":" + drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.join(os.path.relpath(job_args['errorPath'], os.getenv("HOME")))

            jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY + os.sep + os.path.relpath(job_args['workingDirectory'], os.getenv("HOME"))
            jt.nativeSpecification = "-t {time} -p {partition} -A {account} {extra}".format(**job_args)
            if kw.get('email', None):
                jt.email=[kw.get('email')]
            self.app.log.info("Submitting job with native specification {}".format(jt.nativeSpecification))
            self.app.log.info("Working directory: {}".format(jt.workingDirectory))
            self.app.log.info("Output logging: {}".format(jt.outputPath))
            self.app.log.info("Error logging: {}".format(jt.errorPath))
            self._meta.jobid = s.runJob(jt)
            self.app.log.info('Your job has been submitted with id ' + self._meta.jobid)
            if kw.get('saveJobId', False):
                self._save_job_id(**job_args)
            s.deleteJobTemplate(jt)
            s.exit()
        if self.app.pargs.batch:
            self._meta.batch_command.append(command)
        else:
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


def convert_to_drmaa_time(t):
    """Convert time assignment to format understood by drmaa.

    In particular transforms days to hours if provided format is
    d-hh:mm:ss. Also transforms mm:ss to 00:mm:ss.

    :param t: time string

    :returns: converted time string formatted as hh:mm:ss or None if
    time string is malformatted
    """
    if not t:
        return None
    m = re.search("(^[0-9]+\-)?([0-9]+:)?([0-9]+):([0-9]+)", t)
    if not m:
        return None
    days = None
    if m.group(1):
        days = m.group(1).rstrip("-")
    hours = None
    if m.group(2):
        hours = m.group(2).rstrip(":")
    minutes = m.group(3)
    seconds = m.group(4)
    if days:
        hours = 24 * int(days) + int(hours)
    else:
        if not hours:
            hours = "00"
        if len(str(hours)) == 1:
            hours = "0" + hours
        if len(str(minutes)) == 1:
            minutes = "0" + minutes
    t_new = "{}:{}:{}".format(hours, minutes, seconds)
    return t_new

def make_job_template_args(opt_d, **kw):
    """Given a dictionary of arguments, update with kw dict that holds arguments passed to argv.

    :param opt_d: dictionary of option key/value pairs
    :param kw: dictionary of program arguments

    :returns: dictionary of job arguments
    """
    job_args = {}
    job_args['jobname'] = kw.get('jobname', None) or opt_d.get('-J', None) or  opt_d.get('--job-name', None)
    job_args['time'] = kw.get('time', None) or opt_d.get('-t', None) or  opt_d.get('--time', None)
    job_args['time'] = convert_to_drmaa_time(job_args['time'])
    job_args['partition'] = kw.get('partition', None) or opt_d.get('-p', None) or  opt_d.get('--partition', None)
    job_args['account'] = kw.get('account', None) or opt_d.get('-A', None) or  opt_d.get('--account', None)
    job_args['outputPath'] = kw.get('outputPath', None) or opt_d.get('--output', None) or opt_d.get('-o', os.curdir)
    job_args['errorPath'] = kw.get('errorPath', None) or opt_d.get('--error', None) or opt_d.get('-e', os.curdir)
    job_args['workingDirectory'] = kw.get('workingDirectory', None) or opt_d.get('-D', os.curdir) 
    job_args['email'] = kw.get('email', None) or opt_d.get('--mail-user', None) 
    invalid_keys = ["--mail-user", "--mail-type", "-o", "--output", "-D", "--workdir", "-J", "--job-name", "-p", "--partition", "-t", "--time", "-A", "--account", "-e", "--error"]
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
    group.add_argument('--email', help="set user email address", action="store", default=None, type=str)
    group.add_argument('--batch', help="submit jobs as a batch, useful for submitting a number of jobs to the same node", action="store_true", default=False)

def set_distributed_handler(app):
    """
    Overrides the configured command handler if ``--drmaa`` is passed at the
    command line.
    
    :param app: The application object.
    
    """
    if '--drmaa' in app._meta.argv:
        app._meta.cmd_handler = 'distributed'
        app._setup_cmd_handler()

def run_batch_command(app):
    """
    If option 'batch' was set, run commands stored in batch_command
    variable.

    NOTE: currently runs jobs *sequentially*. One could imagine adding
    a pipe to 'gnus parallel' for single-core jobs.

    :param app: The application object.
    """
    if not app.pargs.batch:
        return
    command = ";\n".join(app.cmd._meta.batch_command)
    app.pargs.batch = False
    app.cmd.command([command], **{'platform_args':{}, 'saveJobId':True, 'workingDirectory':os.curdir})

def load():
    """Called by the framework when the extension is 'loaded'."""
    if not os.getenv("DRMAA_LIBRARY_PATH"):
        LOG.debug("No environment variable $DRMAA_LIBRARY_PATH: loading {} failed".format(__name__))
        return
    hook.register('post_setup', add_drmaa_option)
    hook.register('post_setup', add_shared_distributed_options)
    hook.register('pre_run', set_distributed_handler)
    hook.register('post_run', run_batch_command)
    handler.register(DistributedCommandHandler)
