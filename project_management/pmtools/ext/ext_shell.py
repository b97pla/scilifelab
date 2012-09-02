"""Shell extension"""

import subprocess

from cement.core import backend, handler

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
            kwargs = { 'shell': True, 'cwd': cwd}
            if capture:
                kwargs['stderr'] = subprocess.STDOUT
                kwargs['stdout'] = subprocess.PIPE
            p = subprocess.Popen(cmd, **kwargs)
            p_stdout = p.communicate()[0]
            if p.returncode and not ignore_error:
                if capture:
                    self.app.log.error(p_stdout)
                raise Exception("Subprocess return code: %d" % p.returncode)
            if capture:
                return p_stdout
        return self.dry(cmd, runpipe)
            
        
