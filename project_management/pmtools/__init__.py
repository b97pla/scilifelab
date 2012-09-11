"""
Pipeline Management Tools
=========================

Commands
---------

  analysis                
    Manage analysis
  archive                  
    Manage archive
  deliver                  
    Deliver data
  project                  
    Manage projects

"""
__import__('pkg_resources').declare_namespace(__name__)

import os
import sys
import re
import argparse
import textwrap
import subprocess
from cStringIO import StringIO

from cement.core import foundation, controller, handler, backend, output

from pmtools.core import command
from pmtools.core.help import PmHelpFormatter
from pmtools.ext import ext_shell

LOG = backend.minimal_logger(__name__)    

##############################
## Abstract base controller -- for sharing arguments and functions with subclassing controllers
##############################    
class AbstractBaseController(controller.CementBaseController):
    """
    This is an abstract base controller.
    """
    class Meta:
        pass
            
    def _setup(self, base_app):
        self._meta.arguments.append( (['-n', '--dry_run'], dict(help="dry_run - don't actually do anything", action="store_true", default=False)) )
        self._meta.arguments.append((['--force'], dict(help="force execution", action="store_true", default=False)))
        self._meta.arguments.append((['--java_opts'], dict(help="java options", action="store", default="Xmx3g")))
        self._meta.arguments.append((['--input_file'], dict(help="Run on specific input file", default=None)))
        self._meta.arguments.append((['--file_type'], dict(help="file type for globbing", default="")))
        self._meta.arguments.append((['--move'], dict(help="Transfer file with move", default=False, action="store_true")))
        self._meta.arguments.append((['--copy'], dict(help="Transfer file with copy (default)", default=True, action="store_true")))
        super(AbstractBaseController, self)._setup(base_app)

        ## Sometimes read as string, sometimes as list...
        ignore = self.config.get("config", "ignore")
        if type(ignore) == str:
            self.reignore = re.compile(ignore.replace("\n", "|"))
        elif type(ignore) == list:
            self.reignore = re.compile("|".join(ignore))
        self.shared_config = dict()

    def _filtered_ls(self, out):
        """Filter output"""
        def ignore(line):
            return self.reignore.match(line) == None
        return filter(ignore, out)

    def _ls(self, path, filter=False):
        """List contents of path"""
        out = self.app.cmd.command(["ls", path])
        if filter:
            out = self._filtered_ls(out)
        if out:
            self.app._output_data["stdout"].write(out.rstrip())

    def _not_implemented(self, msg=None):
        self.log.warn("FIXME: Not implemented yet")
        if msg != None:
            self.log.warn(msg)
        raise NotImplementedError

    def _obsolete(self, msg):
        self.log.warn("This function is obsolete.")
        self.log.warn(msg)

    def _check_pargs(self, pargs, msg=None):
        """Check that list of pargs are present"""
        for p in pargs:
            if not self.pargs.__getattribute__(p):
                self.log.warn("Required argument '{}' lacking".format(p))
                return False
        return True

    def _assert_config(self, section, label):
        """
        Assert existence of config label. If not present, require
        that the section/label be defined in configuration file.
        """
        if not self.config.has_section(section):
            self.log.error("no such section '{}'; please define in configuration file".format(section))
            sys.exit()
        config_dict = self.config.get_section_dict(section)
        if not config_dict.has_key(label):
            self.log.error("no section '{}' with label '{}' in config file; please define accordingly".format(section, label))
            sys.exit()
        elif config_dict[label] is None:
            self.log.error("config section '{}' with label '{}' set to 'None'; please define accordingly".format(section, label))
            sys.exit()

##############################
## Main controller for all subsubcommands
## Currently does nothing
##############################
class SubController(controller.CementBaseController):
    class Meta:
        pass

    def _setup(self, base_app):
        super(SubController, self)._setup(base_app)

##############################
## Main pm base controller
##############################
class PmController(controller.CementBaseController):
    class Meta:
        label = 'base'
        description = 'Project/pipeline management tools'
        arguments = [
            (['--config'], dict(help="print configuration", action="store_true")),
            (['--config-example'], dict(help="print configuration example", action="store_true")),
            ]

    def _setup(self, app_obj):
        # shortcuts
        super(PmController, self)._setup(app_obj)

    @controller.expose(hide=True)
    def default(self):
        if self.app.pargs.config:
            print "FIXME: show config"
        elif self.app.pargs.config_example:
            print """Configuration example: save as ~/.pm/pm.conf and modify at will.

    [config]
    ignore = slurm*, tmp*

    [archive]
    root = /path/to/archive

    [analysis]
    root = /path/to/illumina

    [log]
    level = INFO
    file = ~/log/pm.log

    [project]
    root = /path/to/projects
    repos = /path/to/repos
        """
        else:
            print self._help_text

##############################
## PmApp
##############################
class PmApp(foundation.CementApp):
    class Meta:
        label = "pm"
        base_controller = PmController
        cmd_handler = ext_shell.ShCommandHandler

    def __init__(self, label=None, **kw):
        super(PmApp, self).__init__(**kw)
        handler.define(command.ICommand)
        self.cmd = None

    def setup(self):
        super(PmApp, self).setup()
        self._setup_cmd_handler()
        self._output_data = dict(stdout=StringIO(), stderr=StringIO())

    def _setup_cmd_handler(self):
        LOG.debug("setting up {}.command handler".format(self._meta.label))
        self.cmd = self._resolve_handler('command', self._meta.cmd_handler)

    def flush(self):
        """Flush output contained in _output_data dictionary"""
        if self._output_data["stdout"].getvalue():
            print self._output_data["stdout"].getvalue()
        if self._output_data["stderr"].getvalue():
            print >> sys.stderr, self._output_data["stderr"].getvalue()
            
##############################
## PmOutputHandler
##############################
class PmOutputHandler(output.CementOutputHandler):
    class Meta:
        label = 'pmout'

    def render(self, data, template = None):
        if data["stdout"].getvalue():
            print >> sys.stdout, data["stdout"].getvalue()
        if data["stderr"].getvalue():
            print >> sys.stderr, data["stderr"].getvalue()
