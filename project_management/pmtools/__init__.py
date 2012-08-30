"""
Pipeline Management Tools
"""
__import__('pkg_resources').declare_namespace(__name__)

import os
import sys
import re
import argparse
import textwrap
import subprocess
from mako.template import Template

from cement.core import foundation, controller, handler, backend
from cement.utils.shell import *


## I personally don't like the default help formatting output from
## argparse which I think is difficult to read
## Currently identical to RawDescriptionHelpFormatter
class PmHelpFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        return ' '.join([indent + line for line in text.splitlines(True)])


##############################
## Abstract base controller -- for sharing arguments and functions with subclassing controllers
## FIXME: make sbatch/drmaa abstract base controller with arguments shared over those functions that need access to slurm (e.g. compress)
##############################    
class AbstractBaseController(controller.CementBaseController):
    """
    This is an abstract base controller.
    """
    class Meta:
        pass
            
    def _setup(self, base_app):
        self._meta.arguments.append( (['-n', '--dry_run'], dict(help="dry_run - don't actually do anything", action="store_true", default=False)) )
        super(AbstractBaseController, self)._setup(base_app)
        self.reignore = re.compile(self.config.get("config", "ignore").replace("\n", "|"))
        self.shared_config = dict()

    def _filtered_ls(self, out):
        """Filter output"""
        def ignore(line):
            return self.reignore.match(line) == None
        return filter(ignore, out)

    def _not_implemented(self, msg=None):
        print "FIXME: Not implemented yet"
        if msg != None:
            print msg

    def _obsolete(self, msg):
        self.log.info("This function is obsolete.")
        self.log.info(msg)

    ## yes or no: http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input
    def query_yes_no(self, question, default="yes"):
        """Ask a yes/no question via raw_input() and return their answer.
        
        "question" is a string that is presented to the user.
        "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).
        
        The "answer" return value is one of "yes" or "no".
        """
        valid = {"yes":True,   "y":True,  "ye":True,
                 "no":False,     "n":False}
        if default == None:
            prompt = " [y/n] "
        elif default == "yes":
            prompt = " [Y/n] "
        elif default == "no":
            prompt = " [y/N] "
        else:
            raise ValueError("invalid default answer: '%s'" % default)

        while True:
            sys.stdout.write(question + prompt)
            choice = raw_input().lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                sys.stdout.write("Please respond with 'yes' or 'no' "\
                                 "(or 'y' or 'n').\n")

    def safe_makedir(self, dname):
        """Make a directory if it doesn't exist"""
        def runpipe():
            if not os.path.exists(dname):
                try:
                    os.makedirs(dname)
                except OSError:
                    if not os.path.isdir(dname):
                        raise
            else:
                self.log.info("Directory %s already exists" % dname)
            return dname
        return self._dry("Make directory %s" % dname, runpipe)

    ## Config helpers - not used?
    def get_dir(self, section, label):
        assert self.config.get(section, label), "no section %s with label %s in config file; please define accordingly" %(section, label)
        d = self.config.get(section,label)
        if not(os.path.exists(d)):
            self.log.warn("no such path %s" % d)
            sys.exit()
        return d

    ## Taken from paver.easy
    ## FIXME: add time stamp (better: make DRY_RUN a log level that only prints to console, for instance by using the interface ILog)
    def _dry(self, message, func, *args, **kw):
        if self.pargs.dry_run:
            print >> sys.stderr, "(DRY_RUN): " + message
            return
        self.log.info(message)
        return func(*args, **kw)

    ## Implement paver-like dry code
    ## NOTE: this is copied from paver.easy, with the slight modification that cmd_args is passed (a list)
    ## FIXME: accept both list (cmd_args) and command (as in paver)? 
    def sh(self, cmd_args, capture=True, ignore_error=False, cwd=None):
        command = " ".join(cmd_args)
        def runpipe():
            kwargs = { 'shell': True, 'cwd': cwd}
            if capture:
                kwargs['stderr'] = subprocess.STDOUT
                kwargs['stdout'] = subprocess.PIPE
            p = subprocess.Popen(command, **kwargs)
            p_stdout = p.communicate()[0]
            if p.returncode and not ignore_error:
                if capture:
                    error(p_stdout)
                raise Exception("Subprocess return code: %d" % p.returncode)
            if capture:
                return p_stdout
        return self._dry(command, runpipe)

    # def sbatch(self, cmd_args, jobname, partition="core"):
    #     command = " ".join(cmd_args)
    #     kw = dict(
    #         mail_type = "ALL",
    #         header = "",
    #         footer = "",
    #         constraint = "",
    #         partition = partition,
    #         jobname = jobname,
    #         mail_user = "per.unneberg@scilifelab.se",
    #         cores = 1,
    #         time = 
    #         )

    #     def runpipe():
    #         pass

    ## FIXME: add sbatch function and templates in case drmaa fails?
    def drmaa(self, cmd_args, jobname, partition="core"):
        if not os.getenv("DRMAA_LIBRARY_PATH"):
            return
        else:
            import drmaa

        if self.pargs.node:
            partition = "node"
        if not self.pargs.uppmax_project:
            self.log.warn("no uppmax id provided; cannot proceed with drmaa command")
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
        return self._dry(command, runpipe)
        ##return self._not_implemented("Implement drmaa code as soon as drmaa library fixed!\nSee AbstractBaseController.drmaa function")
        
    # Copied from cement
    # Modification: - PmHelpFormatter
    #               - only relevant options should be listed for any given command 
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
        self.app.args.formatter_class=PmHelpFormatter

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
            #(['-h','--help'], dict(help="print help", action="store_true")),
            ]

    def _setup(self, app_obj):
        # shortcuts
        super(PmController, self)._setup(app_obj)

    @controller.expose(hide=True)
    def default(self):
        if self.app.pargs.config:
            print "FIXME: show config"
        elif self.app.pargs.config_example:
            print """Configuration example: save as ~/.pm.conf and modify at will.

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
            print __doc__

        




##############################
## PmApp
##############################
class PmApp(foundation.CementApp):
    class Meta:
        label = "pm"
        base_controller = PmController

