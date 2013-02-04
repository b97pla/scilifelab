"""Pm Controller Module"""
import os
import sys
import re
import pprint

from cement.core import interface, handler, controller, backend

from scilifelab.pm.lib.help import PmHelpFormatter
from scilifelab.utils.misc import filtered_output, query_yes_no, filtered_walk

class AbstractBaseController(controller.CementBaseController):
    """
    This is an abstract base controller.

    All controllers should inherit from this class.
    """
    class Meta:
        pattern = ""

    def _setup(self, base_app):
        self._meta.arguments.append( (['-n', '--dry_run'], dict(help="dry_run - don't actually do anything", action="store_true", default=False)))
        self._meta.arguments.append((['--force'], dict(help="force execution", action="store_true", default=False)))
        self._meta.arguments.append((['--verbose'], dict(help="verbose mode", action="store_true", default=False)))
        self._meta.arguments.append((['--java_opts'], dict(help="java options", action="store", default="Xmx3g")))
        super(AbstractBaseController, self)._setup(base_app)
        self.ignore = self.config.get("config", "ignore")
        self.shared_config = dict()

    def _process_args(self):
        pass

    ## Redefine _dispatch so that it also looks for function _process_args
    def _dispatch(self):
        """
        Takes the remaining arguments from self.app.argv and parses for a
        command to dispatch, and if so... dispatches it.
        
        """
        self._add_arguments_to_parser()
        self._parse_args()
        self._process_args()
        
        if not self.command:
            self.app.log.debug("no command to dispatch")
        else:    
            func = self.exposed[self.command]     
            self.app.log.debug("dispatching command: %s.%s" % \
                      (func['controller'], func['label']))

            if func['controller'] == self._meta.label:
                getattr(self, func['label'])()
            else:
                controller = handler.get('controller', func['controller'])()
                controller._setup(self.app)
                getattr(controller, func['label'])()

    def _obsolete(self, msg):
        self.app.log.warn("This function is obsolete.")
        self.app.log.warn(msg)

    def _check_pargs(self, pargs, msg=None):
        """Check that list of pargs are present"""
        for p in pargs:
            if not self.pargs.__getattribute__(p):
                self.app.log.warn("Required argument '{}' lacking".format(p))
                return False
        return True

    def _ls(self, path, filter_output=False):
        """List contents of path"""
        if not os.path.exists(path):
            self.app.log.info("No such path {}".format(path))
            return
        out = self.app.cmd.command(["ls", path])
        if filter_output:
            out = filtered_output(self.ignore, out)
        if out:
            self.app._output_data["stdout"].write(out.rstrip())

class AbstractExtendedBaseController(AbstractBaseController):
    """
    This is an abstract extended base controller.

    All extended controllers should inherit from this class. The main difference to the AbstractBaseController is that this controller adds arguments for compressing and cleaning. 

    """ 
    
    class Meta:
        compress_opt = "-v"
        compress_prog = "gzip"
        compress_suffix = ".gz"
        file_ext = []
        include_dirs = []
        wildcard = []
        root_path = None
        path_id = None

    def _setup(self, base_app):
        super(AbstractExtendedBaseController, self)._setup(base_app)
        group = base_app.args.add_argument_group('Project, sample, and flowcell gorup', 'Options that modify commands to work on specific projects, sample and more.')
        group.add_argument('project', help="Project id", nargs="?", default=None)
        group.add_argument('-f', '--flowcell',help="Flowcell id")
        group.add_argument('-S', '--sample', help="Project sample id. If sample is a file, read file and use sample names within it. Sample names can also be given as full paths to bcbb-config.yaml configuration file.", action="store", default=None, type=str)

        group = base_app.args.add_argument_group('file types', 'Options that set file types.')
        group.add_argument('--sam', help="Workon sam files", default=False, action="store_true")
        group.add_argument('--bam', help="Workon bam files", default=False, action="store_true")
        group.add_argument('--fastq', help="Workon fastq files", default=False, action="store_true")
        group.add_argument('--fastqbam', help="Workon fastq-fastq.bam files", default=False, action="store_true")
        group.add_argument('--pileup', help="Workon pileup files", default=False, action="store_true")
        group.add_argument('--split', help="Workon *-split directories", default=False, action="store_true")
        group.add_argument('--tmp', help="Workon staging (tx) and tmp directories", default=False, action="store_true")
        group.add_argument('--txt', help="Workon txt files", default=False, action="store_true")
        group.add_argument('--glob', help="Workon freetext glob expression. CAUTION: using wildcard expressions will remove *everything* that matches.", default=None, action="store")

        group = base_app.args.add_argument_group('file transfer', 'Options affecting file transfer operations.')
        group.add_argument('--move', help="Transfer file with move", default=False, action="store_true")
        group.add_argument('--copy', help="Transfer file with copy (default)", default=True, action="store_true")
        group.add_argument('--rsync', help="Transfer file with rsync", default=False, action="store_true")

        group = base_app.args.add_argument_group('compression/decompression', 'Options affecting compression/decompression.')
        group.add_argument('--pbzip2', help="Use pbzip2 as compressing device", default=False, action="store_true")
        group.add_argument('--pigz', help="Use pigz as compressing device", default=False, action="store_true")
        group.add_argument('--input_file', help="Run on specific input file", default=None)


    def _process_args(self):
        if self.command in ["compress", "decompress", "clean", "du"]:
            if not self._meta.path_id:
                self.app.log.warn("not running {} on root directory".format(self.command))
                sys.exit()
        if self.pargs.flowcell:
            self._meta.wildcard += [self.pargs.flowcell]
        elif self.command in ["ls"]:
            if not self._meta.path_id:
                self._meta.path_id = ""
        else:
            pass
        ## Setup file patterns to use
        if self.pargs.fastq:
            self._meta.file_ext += [".fastq", "fastq.txt", ".fq"]
        if self.pargs.pileup:
            self._meta.file_ext += [".pileup", "-pileup"]
        if self.pargs.txt:
            self._meta.file_ext += [".txt"]
        if self.pargs.fastqbam:
            self._meta.file_ext += ["fastq-fastq.bam"]
        if self.pargs.sam:
            self._meta.file_ext += [".sam"]
        if self.pargs.bam:
            self._meta.file_ext += [".bam"]
        if self.pargs.split:
            self._meta.file_ext += [".intervals", ".bam", ".bai", ".vcf", ".idx"]
            self._meta.include_dirs += ["realign-split", "variants-split"]
        if self.pargs.tmp:
            self._meta.file_ext += [".idx", ".vcf", ".bai", ".bam", ".idx", ".pdf"]
            self._meta.include_dirs += ["tmp", "tx"]
        if self.pargs.glob:
            self._meta.file_ext += [self.pargs.glob]
            
        ## Setup zip program
        if self.pargs.pbzip2:
            self._meta.compress_prog = "pbzip2"
        elif self.pargs.pigz:
            self._meta.compress_prog = "pigz"

        # Setup transfer options
        if self.pargs.move:
            self.pargs.copy = False
        elif self.pargs.rsync:
            self.pargs.copy = False

        if self._meta.path_id:
            assert os.path.exists(os.path.join(self._meta.root_path, self._meta.path_id)), "no such folder '{}' in {} directory '{}'".format(self._meta.path_id, self._meta.label, self._meta.root_path)
            
    def _filter_fn(self, f):
        if not self._meta.pattern:
            return False
        if self.pargs.flowcell:
            if not re.search(self.pargs.flowcell, f):
                return False
        return re.search(self._meta.pattern, f) != None

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    ## du
    @controller.expose(help="Calculate disk usage")
    def du(self):
        if not self._check_pargs(["project"]):
            return
        out = self.app.cmd.command(["du", "-hs", "{}".format(os.path.join(self._meta.root_path, self._meta.path_id))])
        if out:
            self.app._output_data["stdout"].write(out.rstrip())

    ## clean
    @controller.expose(help="Remove files")
    def clean(self):
        if not self._check_pargs(["project"]):
            return
        self._meta.pattern = "|".join(["{}(.gz|.bz2)?$".format(x) for x in self._meta.file_ext])
        flist = filtered_walk(os.path.join(self._meta.root_path, self._meta.path_id), self._filter_fn, include_dirs=self._meta.include_dirs)
        if len(flist) == 0:
            self.app.log.info("No files matching pattern '{}' found".format(self._meta.pattern))
            return
        if len(flist) > 0 and not query_yes_no("Going to remove {} files ({}...). Are you sure you want to continue?".format(len(flist), ",".join([os.path.basename(x) for x in flist[0:10]])), force=self.pargs.force):
            return
        for f in flist:
            self.app.log.info("removing {}".format(f))
            self.app.cmd.safe_unlink(f)

    def _compress(self, label="compress"):
        if self.pargs.input_file:
            flist = [self.pargs.input_file]
        else:
            flist = filtered_walk(os.path.join(self._meta.root_path, self._meta.path_id), self._filter_fn)

        if len(flist) == 0:
            self.app.log.info("No files matching pattern '{}' found".format(self._meta.pattern))
            return
        if len(flist) > 0 and not query_yes_no("Going to {} {} files ({}...). Are you sure you want to continue?".format(label, len(flist), ",".join([os.path.basename(x) for x in flist[0:10]])), force=self.pargs.force):
            sys.exit()
        for f in flist:
            self.log.info("{}ing {}".format(label, f))
            self.app.cmd.command([self._meta.compress_prog, self._meta.compress_opt, "%s" % f], label, ignore_error=True, **{'workingDirectory':os.path.dirname(f), 'outputPath':os.path.join(os.path.dirname(f), "{}-{}-drmaa.log".format(label, os.path.basename(f)))})

    ## decompress
    @controller.expose(help="Decompress files")
    def decompress(self):
        """Decompress files"""
        if not self._check_pargs(["project"]):
            return
        self._meta.compress_opt = "-dv"
        if self.pargs.pbzip2:
            self._meta.compress_suffix = ".bz2"
        self._meta.pattern = "|".join(["{}{}$".format(x, self._meta.compress_suffix) for x in self._meta.file_ext])
        self._compress(label="decompress")
        
    ## compress
    @controller.expose(help="Compress files")
    def compress(self):
        if not self._check_pargs(["project"]):
            return
        self._meta.compress_opt = "-v"
        self._meta.pattern = "|".join(["{}$".format(x) for x in self._meta.file_ext])
        self._compress()

    ## ls
    @controller.expose(help="List root folder")
    def ls(self):
        if self._meta.path_id == "":
            self._ls(self._meta.root_path, filter_output=True)
        else:
            if self._meta.file_ext:
                pattern = "|".join(["{}$".format(x) for x in self._meta.file_ext])
                flist = filtered_walk(os.path.join(self._meta.root_path, self._meta.path_id), file_filter)
                if flist:
                    self.app._output_data["stdout"].write("\n".join(flist))
            else:
                self._ls(os.path.join(self._meta.root_path, self._meta.path_id))
        
class PmController(controller.CementBaseController):
    """
    Main Pm Controller.

    """
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
            pp = pprint.PrettyPrinter(indent=4)
            out = {k:dict(v) for k,v in self.app.config._sections.iteritems()}
            pp.pprint(out)
        elif self.app.pargs.config_example:
            print """Configuration example: save as ~/.pm/pm.conf and modify at will.

    [config]
    ignore = slurm*, tmp*

    [archive]
    root = /path/to/archive

    [production]
    root = /path/to/production

    [log]
    level = INFO
    file = ~/log/pm.log

    [project]
    root = /path/to/projects
    repos = /path/to/repos
        """
        else:
            print self._help_text
