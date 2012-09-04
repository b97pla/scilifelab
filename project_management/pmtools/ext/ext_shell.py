"""Shell extension"""

import subprocess

from cement.core import backend, handler
from cement.utils import shell

from pmtools.core import command

Log = backend.minimal_logger(__name__)

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

    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        cmd = " ".join(cmd_args)
        def runpipe():
            (stdout, stderr, returncode) = shell.exec_cmd(cmd_args)
            if returncode and not ignore_error:
               if capture:
                   self.app.log.error(stderr)
               raise Exception("Subprocess return code: {}".format(returncode))
            if capture:
               return stdout
            # kwargs = { 'shell': False , 'cwd': cwd}
            # if capture:
            #     kwargs['stderr'] = subprocess.STDOUT
            #     kwargs['stdout'] = subprocess.PIPE
            # p = subprocess.Popen(cmd, **kwargs)
            # p_stdout = p.communicate()[0]
            # if p.returncode and not ignore_error:
            #     if capture:
            #         self.app.log.error(p_stdout)
            #     raise Exception("Subprocess return code: %d" % p.returncode)
            # if capture:
            #     return p_stdout
        return self.dry(cmd, runpipe)
            
        
