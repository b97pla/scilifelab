"""
Pipeline Management Tools

Usage: pm command [options]
"""
import os
import re
import argparse
import textwrap

from cement.core import foundation, controller, handler, backend


## Abstract base controller -- for sharing arguments
class AbstractBaseController(controller.CementBaseController):
    class Meta:
        arguments = [
            (['-n', '--dry_run'], dict(help="dry_run - don't actually do anything")),
            ]

    def _setup(self, base_app):
        super(AbstractBaseController, self)._setup(base_app)
        self.shared_config = dict()

    def _parse_args(self):
        """
        Parse command line arguments and determine a command to dispatch.
        
        """
        # chop off a command argument if it matches an exposed command
        if len(self.app.argv) > 0 and not self.app.argv[0].startswith('-'):
            # translate dashes back to underscores
            cmd = re.sub('-', '_', self.app.argv[0])
            if cmd in self.exposed:
                visible = {cmd : self.visible[cmd]}
                self.visible = visible
                self.command = cmd
                self.app.argv.pop(0)
            else:
                for label in self.exposed:
                    func = self.exposed[label]
                    if self.app.argv[0] in func['aliases']:
                        self.command = func['label']
                        self.app.argv.pop(0)
                        break
                        
        self.app.args.description = self._help_text
        self.app.args.usage = self._usage_text
        self.app.args.formatter_class=argparse.RawDescriptionHelpFormatter

        self.app._parse_args()
        self.pargs = self.app.pargs

    @property
    def _usage_text(self):
        """
        Returns the usage text displayed when '--help' is passed.
        
        """
        if self == self.app._meta.base_controller:
            txt = "%s <CMD> -opt1 --opt2=VAL [arg1] [arg2] ..." % \
                self.app.args.prog
        elif self.command != "default":
            txt = "%s %s %s -opt1 --opt2=VAL [arg1] [arg2] ..." % \
                (self.app.args.prog, self._meta.label, self.command)
        else:
            txt = "%s %s <CMD> -opt1 --opt2=VAL [arg1] [arg2] ..." % \
                  (self.app.args.prog, self._meta.label)
        return txt

    @property
    def _help_text(self):
        """
        Returns the help text displayed when '--help' is passed.
        
        """
        cmd_txt = ''
        
        # hack it up to keep commands in alphabetical order
        sorted_labels = []
        
        for label in list(self.visible.keys()):
            old_label = label
            label = re.sub('_', '-', label)
            sorted_labels.append(label)
            
            if label != old_label:
                self.visible[label] = self.visible[old_label]
                del self.visible[old_label]
        sorted_labels.sort()

        for label in sorted_labels:
            func = self.visible[label]
            if len(func['aliases']) > 0:
                cmd_txt = cmd_txt + "  %s (aliases: %s)\n" % \
                            (label, ', '.join(func['aliases']))
            else:
                cmd_txt = cmd_txt + "  %s\n" % label
            
            if func['help']:
                cmd_txt = cmd_txt + "    %s\n\n" % func['help']
            else:
                cmd_txt = cmd_txt + "\n"
        if len(cmd_txt) > 0:
            if self.command == "default":
                txt = '''%s
                
commands:
                
%s


''' % (self._meta.description, cmd_txt)
            else:
                txt = '''
%s
''' % (cmd_txt)
        else:
            txt = self._meta.description
        return textwrap.dedent(txt)        


## Main controller for all subsubcommands
## Currently does nothing
class SubSubController(controller.CementBaseController):
    class Meta:
        pass

    def _setup(self, base_app):
        super(SubSubController, self)._setup(base_app)


##############################
## Main pm base controller
##############################
class PmController(controller.CementBaseController):
    class Meta:
        label = 'pm'
        description = 'Project/pipeline management tools'
        config_defaults = dict(
                 config="%s/.pm/pm.conf" % os.getenv("HOME")
            )
        arguments = [
            (['--config'], dict(help = "pm config file [default ~/.pm/pm.conf]", action="store", dest="config"))
            ]


    @controller.expose(hide=True)
    def default(self):
        defaults =  backend.defaults('pm')
        print defaults['pm']

    def _parse_args(self):
        """
        Parse command line arguments and determine a command to dispatch.
        
        """
        if len(self.app.argv) == 0:
            self.app.argv.append("-h")

        # chop off a command argument if it matches an exposed command
        if len(self.app.argv) > 0 and not self.app.argv[0].startswith('-'):
            
            # translate dashes back to underscores
            cmd = re.sub('-', '_', self.app.argv[0])
            if cmd in self.exposed:
                self.command = cmd
                self.app.argv.pop(0)
            else:
                for label in self.exposed:
                    func = self.exposed[label]
                    if self.app.argv[0] in func['aliases']:
                        self.command = func['label']
                        self.app.argv.pop(0)
                        break
        self.app.args.description = self._help_text
        self.app.args.usage = self._usage_text
        self.app.args.formatter_class=argparse.RawDescriptionHelpFormatter

        self.app._parse_args()
        self.pargs = self.app.pargs
