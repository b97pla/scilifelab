"""Pm production module"""

import os
import re
import yaml
import platform
import time
import glob
import shutil
import subprocess

from datetime import datetime

from cement.core import controller
from scilifelab.pm.core.controller import AbstractExtendedBaseController
from scilifelab.utils.misc import query_yes_no, filtered_walk, last_lines
from scilifelab.bcbio import prune_pp_platform_args
from scilifelab.bcbio.flowcell import Flowcell
from scilifelab.bcbio.status import status_query
from scilifelab.utils.string import strip_extensions
from scilifelab.pm.core.bcbio import BcbioRunController
from scilifelab.utils.timestamp import utc_time
from scilifelab.db.statusdb import FlowcellRunMetricsConnection

FINISHED_FILE = "FINISHED_AND_DELIVERED"
REMOVED_FILE = "FINISHED_AND_REMOVED"

## Main production controller
class ProductionController(AbstractExtendedBaseController, BcbioRunController):
    """
    Functionality for production management.
    """
    class Meta:
        label = 'production'
        description = 'Manage production'

    def _setup(self, base_app):
        ## Adding arguments to existing groups requires setting up parent class first
        super(ProductionController, self)._setup(base_app)
        group = [x for x in self.app.args._action_groups if x.title == 'file transfer'][0]
        group.add_argument('--from_pre_casava', help="Transfer file with move", default=False, action="store_true")
        group.add_argument('--to_pre_casava', help="Use pre-casava directory structure for delivery", action="store_true", default=False)
        group.add_argument('--transfer_dir', help="Transfer data to transfer_dir instead of sample_prj dir", action="store", default=None)
        base_app.args.add_argument('--brief', help="Output brief information from status queries", action="store_true", default=False)

    def _process_args(self):
        # Set root path for parent class
        ## FIXME: use abspath?
        self._meta.root_path = self.app.config.get("production", "root")
        assert os.path.exists(self._meta.root_path), "No such directory {}; check your production config".format(self._meta.root_path)
        ## Set path_id for parent class
        if self.pargs.flowcell:
            self._meta.path_id = self.pargs.flowcell
        if self.pargs.project:
            self._meta.path_id = self.pargs.project
        ## Temporary fix for pre-casava directories
        if self.pargs.from_pre_casava:
            self._meta.path_id = self.pargs.flowcell
        ## This is a bug; how will this work when processing casava-folders?!?
        ## I need to set this so as not to upset productioncontrollers process_args
        if self.command == "hs_metrics":
            self._meta.path_id = self.pargs.flowcell if self.pargs.flowcell else self.pargs.project
        super(ProductionController, self)._process_args()

    @controller.expose(hide=True)
    def default(self):
        print self._help_text

    @controller.expose(help="Query the status of flowcells, projects, samples"\
                           " that are organized according to the CASAVA file structure")
    def status_query(self):
        if not self._check_pargs(["project", "flowcell"]):
            return
        status_query(self.app.config.get("archive", "root"), self.app.config.get("production", "root"), self.pargs.flowcell, self.pargs.project, brief=self.pargs.brief)

    def _from_casava_structure(self):
        """Get information from casava structure"""
        if not self._check_pargs(["project"]):
            return
        fc_list = []
        pattern = "-bcbb-config.yaml$"
        def bcbb_yaml_filter(f):
            return re.search(pattern, f) != None
        samples = filtered_walk(os.path.join(self._meta.root_path, self._meta.path_id), bcbb_yaml_filter)
        for s in samples:
            fc = Flowcell(s)
            fc_new = fc.subset("sample_prj", self.pargs.project)
            fc_new.collect_files(os.path.dirname(s))
            fc_list.append(fc_new)
        return fc_list

    def _to_casava_structure(self, fc):
        transfer_status = {}
        outdir_pfx = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project, "data"))
        if self.pargs.transfer_dir:
           outdir_pfx = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.transfer_dir, "data"))
        for sample in fc:
            key = "{}_{}".format(sample['lane'], sample['sequence'])
            sources = {"files": self._prune_sequence_files(sample['files']), "results":sample['results']}
            outdir = os.path.join(outdir_pfx, sample['name'], fc.fc_id())
            dirs = {"data":os.path.abspath(os.path.join(outdir_pfx, sample['name'], fc.fc_id())),
                    "intermediate":os.path.abspath(os.path.join(outdir_pfx, sample['name'], fc.fc_id()))}
            self._make_output_dirs(dirs)
            fc_new = fc.subset("lane", sample['lane']).subset("name", sample['name'])
            targets = {"files": [src.replace(fc.path, dirs["data"]) for src in sources['files']],
                       "results": [src.replace(fc.path, dirs["intermediate"]) for src in sources['results']]}

            fc_new.lane_files = dict((k, [os.path.join(outdir, os.path.basename(x)) for x in v]) for k,v in fc_new.lane_files.items())
            fc_new.set_entry(key, 'files', targets['files'])
            fc_new.set_entry(key, 'results', targets['results'])
            ## Copy sample files - currently not doing lane files
            self._transfer_files(sources, targets)
            self.app.cmd.write(os.path.join(dirs["data"], "{}-bcbb-pm-config.yaml".format(sample['name'])), fc_new.as_yaml())
            transfer_status[sample['name']] = {'files':len(sources['files']), 'results':len(sources['results'])}
        ## Rewrite platform_args; only keep time, workdir, account, partition, outpath and jobname
        pattern = "-post_process.yaml$"
        def pp_yaml_filter(f):
            return re.search(pattern, f) != None
        ppfiles = filtered_walk(dirs["data"], pp_yaml_filter)
        for pp in ppfiles:
            self.app.log.debug("Rewriting platform args for {}".format(pp))
            with open(pp, "r") as fh:
                conf = yaml.load(fh)
            if not conf:
                self.app.log.warn("No configuration for {}".format(pp))
                continue
            newconf = prune_pp_platform_args(conf)
            if newconf == conf:
                continue
            self.app.cmd.safe_unlink(pp)
            self.app.cmd.write(pp, yaml.safe_dump(newconf, default_flow_style=False, allow_unicode=True, width=1000))

        # Write transfer summary
        self.app._output_data["stderr"].write("Transfer summary\n")
        self.app._output_data["stderr"].write("{:<18}{:>18}{:>18}\n".format("Sample","Transferred files", "Results"))
        for k, v in transfer_status.iteritems():
            self.app._output_data["stderr"].write("{:<18}{:>18}{:>18}\n".format(k, v['files'], v['results']))

    def _to_pre_casava_structure(self, fc):
        dirs = {"data":os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project, "data", fc.fc_id())),
                "intermediate":os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.project, "intermediate", fc.fc_id()))}
        if self.pargs.transfer_dir:
           dirs["data"] = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.transfer_dir, "data", fc.fc_id()))
           dirs["intermediate"] = os.path.abspath(os.path.join(self.app.config.get("project", "root"), self.pargs.transfer_dir, "intermediate", fc.fc_id()))
        self._make_output_dirs(dirs)
        fc_new = fc
        for sample in fc:
            key = "{}_{}".format(sample['lane'], sample['sequence'])
            sources = {"files": self._prune_sequence_files(sample['files']), "results":sample['results']}
            targets = {"files": [src.replace(fc.path, dirs["data"]) for src in sources['files']],
                       "results": [src.replace(fc.path, dirs["intermediate"]) for src in sources['results']]}
            fc_new.set_entry(key, 'files', targets['files'])
            fc_new.set_entry(key, 'results', targets['results'])
            ## FIX ME: lane file gathering
            ## fc_new.lane_files = dict((k,[x.replace(indir, outdir) for x in v]) for k,v in fc_new.lane_files.items())
            ## Copy sample files - currently not doing lane files
            self._transfer_files(sources, targets)
        self.app.cmd.write(os.path.join(dirs["data"], "project_run_info.yaml"), fc_new.as_yaml())

    def _from_pre_casava_structure(self):
        if not self._check_pargs(["project", "flowcell"]):
            return
        fc = Flowcell()
        fc.load([os.path.join(x, self.pargs.flowcell) for x in [self.config.get("archive", "root"), self.config.get("production", "root")]])
        indir = os.path.join(self.config.get("production", "root"), self.pargs.flowcell)
        if not fc:
            self.log.warn("No run information available for {}".format(self.pargs.flowcell))
            return
        fc_new = fc.subset("sample_prj", self.pargs.project)
        fc_new.collect_files(indir)
        return fc_new

    def _prune_sequence_files(self, flist):
        """Sometimes fastq and fastq.gz files are present for the same
        read. Make sure only one file is used, giving precedence to
        the zipped file.

        :param flist: list of files

        :returns: pruned list of files
        """
        samples = {k[0]:{} for k in [strip_extensions(x, ['.gz']) for x in flist]}
        tmp = {samples[str(k[0])].update({str(k[1]):True}) for k in [strip_extensions(x, ['.gz']) for x in flist]}
        return ["{}.gz".format(k) if v.get('.gz', None) else k for k, v in samples.items()]

    def _make_output_dirs(self, dirs):
        if not os.path.exists(dirs["data"]):
            self.app.cmd.safe_makedir(dirs["data"])
        if not os.path.exists(dirs["intermediate"]):
            self.app.cmd.safe_makedir(dirs["intermediate"])

    def _transfer_files(self, sources, targets):
        for src, tgt in zip(sources['files'] + sources['results'], targets['files'] + targets['results']):
            if not os.path.exists(os.path.dirname(tgt)):
                self.app.cmd.safe_makedir(os.path.dirname(tgt))
            self.app.cmd.transfer_file(src, tgt)

    @controller.expose(help="Transfer data")
    def transfer(self):
        if not self.pargs.from_pre_casava and self.pargs.to_pre_casava:
            self.app.log.warn("not delivering from casava input to pre_casava output")
            return
        ## Collect files depending on input structure
        if self.pargs.from_pre_casava:
            fc = self._from_pre_casava_structure()
        else:
            fc = self._from_casava_structure()
        if not fc:
            return
        ## Organize output file names depending on output structure
        if self.pargs.to_pre_casava:
            self._to_pre_casava_structure(fc)
        else:
            if isinstance(fc, list):
                for f in fc:
                    self._to_casava_structure(f)
            else:
                self._to_casava_structure(fc)

    ## Command for touching file that indicates finished samples
    @controller.expose(help="Touch finished samples. Creates a file FINISHED_AND_DELIVERED with a utc time stamp.")
    def touch_finished(self):
        if not self._check_pargs(["project", "sample"]):
            return
        if os.path.exists(self.pargs.sample) and os.path.isfile(self.pargs.sample):
            with open(self.pargs.sample) as fh:
                slist = [x.rstrip() for x in fh.readlines()]
        else:
            slist = [self.pargs.sample]
        for s in slist:
            spath = os.path.join(self._meta.root_path, self._meta.path_id, s)
            if not os.path.exists(spath):
                self.app.log.warn("No such path {}; skipping".format(spath))
                continue
            rsync_src = os.path.join(self._meta.root_path, self._meta.path_id, s) + os.sep
            rsync_tgt = os.path.join(self.app.config.get("runqc", "root"), self.pargs.project, s) + os.sep
            cl = ["rsync {} {} {}".format(self.app.config.get("runqc", "rsync_sample_opts"), rsync_src, rsync_tgt)]
            self.app.log.info("Checking if runqc uptodate with command '{}'".format(" ".join(cl)))
            out = self.app.cmd.command(cl, **{'shell':True})
            if not self.pargs.dry_run and not out.find("total size is 0"):
                self.app.log.info("Some files need to be updated. Rsync output:")
                print "********"
                print out
                print "********"
                continue
            if not query_yes_no("Going to touch file {} for sample {}; continue?".format(FINISHED_FILE, s), force=self.pargs.force):
                continue
            self.app.log.info("Touching file {} for sample {}".format(FINISHED_FILE, s))
            with open(os.path.join(spath, FINISHED_FILE), "w") as fh:
                t_utc = utc_time()
                fh.write(t_utc)

    ## Command for removing samples that have a FINISHED_FILE flag
    @controller.expose(help="Remove finished samples for a project. Searches for FINISHED_AND_DELIVERED and removes sample contents if file is present.")
    def remove_finished(self):
        if not self._check_pargs(["project"]):
            return
        # Don't filter out files
        def filter_fn(f):
            return True
        slist = os.listdir(os.path.join(self._meta.root_path, self._meta.path_id))
        for s in slist:
            spath = os.path.join(self._meta.root_path, self._meta.path_id, s)
            if not os.path.isdir(spath):
                continue
            if not os.path.exists(os.path.join(spath, FINISHED_FILE)):
                self.app.log.info("Sample {} not finished; skipping".format(s))
                continue
            flist = filtered_walk(spath, filter_fn)
            dlist = filtered_walk(spath, filter_fn, get_dirs=True)
            if os.path.exists(os.path.join(spath, REMOVED_FILE)):
                self.app.log.info("Sample {} already removed; skipping".format(s))
                continue
            if len(flist) > 0 and not query_yes_no("Will remove directory {} containing {} files; continue?".format(s, len(flist)), force=self.pargs.force):
                continue
            self.app.log.info("Removing {} files from {}".format(len(flist), spath))
            for f in flist:
                if f == os.path.join(spath, FINISHED_FILE):
                    continue
                self.app.cmd.safe_unlink(f)
            self.app.log.info("Removing {} directories from {}".format(len(dlist), spath))
            for d in sorted(dlist, reverse=True):
                self.app.cmd.safe_rmdir(d)
            if not self.pargs.dry_run:
                with open(os.path.join(spath, REMOVED_FILE), "w") as fh:
                    t_utc = utc_time()
                    fh.write(t_utc)

    @controller.expose(help="Cleans up the storage systems for production environment. " \
            "It will distinguish between primary storage systems (i.e production NAS) " \
            "and analysis machines (i.e b5/UPPMAX)")
    def storage_cleanup(self):
        storage_conf = self.app.config.get_section_dict('storage')
        db_info = self.app.config.get_section_dict('db')
        f_conn = FlowcellRunMetricsConnection(username=db_info.get('user'),
                                              password=db_info.get('password'),
                                              url=db_info.get('url'))
        servers = [server for server in storage_conf.keys()]
        server = platform.node().split('.')[0].lower()
        if server in servers:
            self.app.log.info("Performing cleanup on production server \"{}\"...".format(server))
            dirs = [d.lstrip() for d in storage_conf.get(server).split(',')]

            #Collect old runs (> 30 days in nosync folder) to remove
            old_runs = []
            for d in dirs:
                nosync_dir = os.path.join(d, 'nosync')
                for fc in glob.glob(os.path.join(nosync_dir, '1*')):
                    fc_name = os.path.basename(fc)
                    #Check that there is no check file indicating to not remove the run
                    if not os.path.exists(os.path.join(fc, 'no_remove.txt')):
                        stats = os.stat(os.path.join(fc, 'RTAComplete.txt'))
                        mod_time = datetime.now() - datetime.fromtimestamp(stats.st_mtime)
                        if mod_time.days >= 30:
                            old_runs.append(fc)
                    else:
                        self.app.log.warn("no_remove.txt file found in {}, skipping run".format(fc_name))

            #NAS servers
            if 'nas' in server:
                #Collect newly finished runs
                fc_list = []
                for d in dirs:
                    for fc in glob.glob(os.path.join(d, '1*')):
                        if os.path.exists(os.path.join(fc, 'RTAComplete.txt')):
                            fc_list.append(fc)

                #Move to nosync
                retries = 5
                for fc in fc_list:
                    fc_name = os.path.basename(fc)
                    while retries:
                        if 'Finished' in last_lines(storage_conf.get('lsyncd_log'), 1)[0]:
                            break
                        retries -= 1
                        time.sleep(3)
                    if retries:
                        self.app.log.info("lsyncd process seems to be up to speed, and run {} " \
                                "is finished, moving it to nosync".format(fc_name))
                        shutil.move(fc, os.path.join(os.path.dirname(fc), 'nosync'))
                        #Touch RTAComplete.txt file to that the modification date is the date when
                        #it was moved to nosync
                        open(os.path.join(os.path.dirname(fc), 'nosync', os.path.basename(fc), 'RTAComplete.txt'), 'w').close()
                        fc_db_id = f_conn.id_view.get(fc_name)
                        f_conn.set_storage_status(fc_db_id, 'NAS_nosync')
                    else:
                        self.app.log.warn("lsyncd process doesn't seem to be finished. " \
                                "Skipping run {}".format(os.path.basename(fc)))

                #Remove old runs
                for fc in old_runs:
                    fc_name = os.path.basename(fc)
                    #Check that the run has been archived in swestore before removing permanently
                    if fc_name in f_conn.get_storage_status('swestore_archived').keys():
                        self.app.log.info("Run {} has been in nosync for more than 30 days " \
                            "and is archived in swestore. Permanently removing it from the NAS".format(fc_name))
                        shutil.rmtree(fc)
                    else:
                        self.app.log.warn("Run {} has been in nosync for more than 30 " \
                            "days, but has not yet been archived in swestore. " \
                            "Not removing, please check it".format(fc_name))

            #Processing servers (b5)
            else:
                #Collect finished runs
                fc_list = []
                for d in dirs:
                    for fc in glob.glob(os.path.join(d, '1*')):
                        if os.path.exists(os.path.join(fc, 'second_read_processing_completed.txt')):
                            fc_list.append(fc)

                #Move to nosync
                for fc in fc_list:
                    fc_name = os.path.basename(fc)
                    self.app.log.info("Moving run {} to nosync".format(fc_name))
                    shutil.move(fc, os.path.join(os.path.dirname(fc), 'nosync'))

                #Remove old runs
                for fc in old_runs:
                    fc_name = os.path.basename(fc)
                    self.app.log.info("Run {} has been in nosync for more than 30 " \
                        "days, permanently removing it from {}".format(fc_name, server))
                    shutil.rmtree(fc)
        else:
            self.app.log.warn("You're running the cleanup functionality in {}. But this " \
                    "server doen't seem to be on your pm.conf file. Are you on the correct server?".format(server))

    @controller.expose(help="Syncs a run to UPPMAX swestre_staging dir. It uses a regular " \
            "rsync command without preserving modification times.")
    def sync_run(self):
        storage_conf = self.app.config.get_section_dict('storage')
        archive_conf = self.app.config.get_section_dict('archive')
        swestore_dir = self.app.config.get_section_dict('archive').get('swestore_staging')
        servers = [server for server in storage_conf.keys()]
        server = platform.node().split('.')[0].lower()
        flowcell = self.pargs.flowcell
        if not flowcell:
            self.app.log.error("Flowcell parameter is required.")
            return
        if server in servers:
            #Find the run of the run directories
            run_dir = ''
            for d in [d.lstrip() for d in storage_conf.get(server).split(',')]:
                if os.path.exists(os.path.join(d, flowcell)):
                    run_dir = os.path.join(d, flowcell)
            if not run_dir:
                self.app.log.error("Run {} not found on the server. ")
                return

            cl = ['rsync',
                  '--recursive',
                  run_dir,
                  '{}@{}:{}'.format(archive_conf.get('user'),
                                    archive_conf.get('server'),
                                    archive_conf.get('swestore_staging'))]
            try:
                subprocess.check_call(cl)
            except subprocess.CalledProcessError():
                self.app.log.error("Rsync of run {} failed, please check that you have " \
                        "the correct username and server in your pm.conf file, " \
                        "current ones are {} and {}".format(flowcell, archive_conf.get('user'),
                                                            archive_conf.get('server')))
                return

        else:
            self.app.log.warn("You're running the cleanup functionality in {}. But this " \
                    "server doen't seem to be on your pm.conf file. Are you on the correct server?".format(server))
