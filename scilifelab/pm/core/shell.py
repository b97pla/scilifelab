"""Shell core module"""
import os
import psutil
import subprocess

from cement.core import backend, handler
from cement.utils import shell

from scilifelab.pm.core import command

LOG = backend.minimal_logger(__name__)

class ShCommandHandler(command.CommandHandler):
    """ 
    This class is an implementation of the :ref:`ICommand
    <pmtools.core.command>` interface.
    """
    class Meta:
        """Handler meta-data"""
        
        interface = command.ICommand
        """The interface that this class implements."""

        label = 'shell'
        """The string identifier of this handler."""

        pid = None
        """The process id"""

    def monitor(self, work_dir, idfile="PID"):
        """Check for existing process"""
        return self._monitor_process(work_dir, idfile="PID")

    def _save_process_id(self, work_dir, idfile="PID"):
        """Save process to file"""
        PIDFILE = os.path.join(work_dir, idfile)
        with open(PIDFILE, "w") as fh:
            fh.write(pid)
        
    def _monitor_process(self, work_dir, idfile="PID"):
        """Monitor existing process"""
        PIDFILE = os.path.join(work_dir, idfile)
        if not os.path.exists(PIDFILE):
            return
        with open(PIDFILE) as fh:
            pid = fh.read()
        if psutil.pid_exists(pid):
            self.app.log.warn("pid {} still running; please terminate job before proceeding".format(pid))
            return True
        else:
            return False
            

    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        cmd = " ".join(cmd_args)
        def runpipe():
            (stdout, stderr, returncode) = shell.exec_cmd(cmd, kw.get('shell', False))
            if returncode and not ignore_error:
               if capture:
                   self.app.log.error(stderr)
               raise Exception("Subprocess return code: {}".format(returncode))
            if capture:
               return stdout
        return self.dry(cmd, runpipe)
            
        
