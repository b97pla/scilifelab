#!/usr/bin/env python
"""Data management based on project identifiers

Usage:
  project_management.py     <YAML config file> <flowcell_dir> <project_dir>
                            [<YAML run information>
                             --flowcell_alias=<flowcell_alias> --project_desc=<project_desc>
                             --install_data --move_data --move_and_relink --symlink
                             --only_install_run_info --customer_delivery --barcode_id_to_name
                             --no_full_names --dry_run --verbose]


Given a directory with demultiplexed flow cell data and a project id,
project specific files will be copied to a project directory. The
optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. See
automated_initial_analysis.py for details.

For a multiproject run_info file, only a subset of the lanes can be
used. The run_info file is therefore pruned, and the pruned file project_run_info.yaml 
is output to the project directory. The pruning is based on the options
<project_desc>. 

In order to avoid downstream demultiplexing, a file sample_project_run_info.yaml 
is created, in which each sample is put in a separate lane.

Options:
  -a, --flowcell_alias=<flowcell alias>         By default, samples are moved to a directory named
                                                <flowcell_dir>. This option links the output directory to
                                                <flowcell_alias>.
  -y, --project_desc=<project_desc>             Project description in sample_prj field of run_info file, or ALL.
  -c, --customer_delivery                       Deliver data to customer. Delivers all demuxed fastq data and results
                                                to one directory <project_desc/flowcell_dir>.
  -b, --barcode_id_to_name                      Convert barcode ids to sample names
  -i, --only_install_run_info                   Only install pruned run_info file.
  -f, --no_full_names                           Don't use full barcode name, including "_index"
  -m, --move_data                               Move data instead of copying
  -M, --move_and_relink                         Move data from source to target, and link from target to source
  -l, --symlink                                 Link data instead of copying
  -n, --dry_run                                 Don't do anything samples, just list what will happen
  -v, --verbose                                 Print some more information
"""

import os
import sys
import re
from optparse import OptionParser

import yaml
import glob
import shutil
import logbook
from itertools import izip

from bcbio.log import logger, setup_logging
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.lane import get_flowcell_id
from bcbio.pipeline.fastq import get_single_fastq_files#, get_fastq_files
from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.flowcell import Flowcell, Lane, get_sample_name
from bcbio import utils
from bcbio.pipeline.config_loader import load_config

class PostProcessedFlowcell(Flowcell):
    """A class for managing information about a post processed flowcell.

    Adds functionality for keeping track of where analysis results
    are stored, and adds information about fc id and fc alias for data delivery"""

    def __init__(self, fc_name, fc_date, data, fc_dir=None, fc_alias=None, fc_results_dir=None):
        self.set_fc_dir(fc_dir)
        self.set_fc_date(fc_date)
        self.set_fc_name(fc_name)
        self.set_lanes(data)
        if fc_alias is None:
            fc_alias = self.get_fc_id()
        self.set_fc_alias(fc_alias)
        if fc_results_dir is None:
            fc_results_dir = fc_dir
        self.set_fc_results_dir(fc_results_dir)

    def get_files(self):
        fl = []
        for l in self.get_lanes():
            fl.append(l.get_files())
        return fl

    def get_fc_id(self):
        return "%s_%s" % (self.get_fc_date(), self.get_fc_name())

    def get_fc_results_dir(self):
        return self.fc_results_dir
    def set_fc_results_dir(self,fc_results_dir):
        self.fc_results_dir = fc_results_dir

    def get_fc_alias(self):
        return self.fc_alias
    def set_fc_alias(self,fc_alias):
        self.fc_alias = fc_alias

    def prune_to_project(self,project, exclude_unmatched=True):
        """Return a new PostProcessedFlowcell object just containing the lanes and samples belonging to a specific project"""
        lanes = []
        fc = None
        for lane in self.get_lanes():
            l = lane.prune_to_project(project, exclude_unmatched)
            if (l):
                lanes.append(l.to_structure())
        if (len(lanes)):
            fc = PostProcessedFlowcell(self.get_fc_name(),self.get_fc_date(),lanes,self.get_fc_dir())
        return fc

    def __str__(self):
        s =  """PostProcessedFlowcell: %s
    fc_dir:          %s
    fc_date:         %s
    fc_name:         %s
    fc_id:           %s
    fc_alias:        %s
    fc_results_dir:  %s

    Number of lanes: %s
""" % (self.get_fc_id(), self.get_fc_dir(), self.get_fc_date(), self.get_fc_name(), self.get_fc_id(), self.get_fc_alias(), self.get_fc_results_dir(), len(self.lanes))
        return s

def main(config_file, fc_dir, project_dir, run_info_yaml=None, fc_alias=None, project_desc=None, lanes=None, barcodes=None):
    config = load_config(config_file)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(project_dir, "log")
    setup_logging(config)

    if project_desc is None and lanes is None:
        logger.error("No project description or lanes provided: cannot deliver files without this information")
        sys.exit()

    if options.customer_delivery and not fc_alias == "":
        logger.info("INFO: Ignoring flowcell_alias when doing customer_delivery")

    fc_dir = os.path.abspath(fc_dir)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fp = open(run_info_yaml)
    run_info_structure = yaml.load(fp)
    original_fc = PostProcessedFlowcell(fc_name, fc_date, run_info_structure, fc_dir=fc_dir, fc_results_dir=fc_dir)
    pruned_fc = original_fc.prune_to_project(project_desc, exclude_unmatched=True)
    if pruned_fc is None or len(pruned_fc.get_lanes()) == 0:
        if not project_desc is None:
            logger.error("No lanes found with matching description %s: please check your flowcell run information" % project_desc)
            print >> sys.stderr, "Available projects: \n\t%s" %  ("\n\t".join(original_fc.get_project_names()))
            sys.exit()
        if not lanes  is None:
            logger.error("No lanes found with numbers %s: please check your flowcell run information" % " ".join(lanes))
            sys.exit()
    # Set up a raw data flowcell that contains the delivery information for raw data (demuxed fastq data)
    rawdata_fc = PostProcessedFlowcell(pruned_fc.get_fc_name(), pruned_fc.get_fc_date(), pruned_fc.to_structure()['details'], fc_alias=fc_alias)
    rawdata_fc.set_fc_dir(os.path.abspath(os.path.join(project_dir, "nobackup/data", rawdata_fc.get_fc_id())))
    analysis_fc = PostProcessedFlowcell(pruned_fc.get_fc_name(), pruned_fc.get_fc_date(), pruned_fc.to_structure()['details'], fc_alias=fc_alias)
    analysis_fc.set_fc_dir(os.path.abspath(os.path.join(project_dir, "nobackup/intermediate", rawdata_fc.get_fc_id())))

    # If customer delivery setup some special options
    if options.customer_delivery:
        rawdata_fc.set_fc_dir(os.path.abspath(os.path.join(project_dir, project_desc, rawdata_fc.get_fc_id())))
        rawdata_fc.set_fc_alias(rawdata_fc.get_fc_id())
        analysis_fc = rawdata_fc
    _make_delivery_directory(rawdata_fc)
    _make_delivery_directory(analysis_fc)
    run_main(pruned_fc, rawdata_fc, analysis_fc)

def run_main(pruned_fc, rawdata_fc, analysis_fc):
    rawdata_fc.set_lanes([])
    for lane in pruned_fc.get_lanes():
        new_lane = process_lane(lane, pruned_fc, rawdata_fc, analysis_fc)
        rawdata_fc.add_lane(new_lane)
    _save_run_info(rawdata_fc, "project_run_info.yaml")
    _sample_based_run_info(rawdata_fc)

def process_lane(lane, pruned_fc, rawdata_fc, analysis_fc):
    """Models bcbio process lane"""
    multiplex = lane.get_samples()
    logger.info("Processing lane %s; reference genome %s" %
             (lane.get_name(), lane.get_genome_build()))
    if multiplex:
        logger.debug("Project %s is multiplexed as: %s" % (lane.get_name(), multiplex))
    fq = _get_barcoded_fastq_files(lane, multiplex, pruned_fc.get_fc_date(), pruned_fc.get_fc_name(), pruned_fc.get_fc_dir())

    ## Move data along with fastq files
    fc_data_dir = rawdata_fc.get_fc_dir()
    _make_dir(fc_data_dir, "data delivery directory")
    if options.install_data:
        data, fastqc = _get_analysis_results(pruned_fc, lane)
        _deliver_data(data, fastqc, analysis_fc.get_fc_dir())
    fastq_targets = list()
    for fqpair in fq:
        for fastq_src in fqpair:
            fastq_tgt = fastq_src
            if options.customer_delivery or options.barcode_id_to_name:
                fastq_tgt = _convert_barcode_id_to_name(multiplex, rawdata_fc.get_fc_name(), fastq_src)
                fastq_tgt.replace("_nophix_", "_")
            _deliver_fastq_file(fastq_src, os.path.basename(fastq_tgt), fc_data_dir)
            fastq_targets.append(os.path.join(fc_data_dir, os.path.basename(fastq_tgt)))
    lane.set_files(fastq_targets)
    return lane

def _get_barcoded_fastq_files(lane, multiplex, fc_date, fc_name, fc_dir=None):
    fq = list()
    bc_dir = "%s_%s_%s_nophix_barcode" % (lane.get_name(), fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    if multiplex is None:
        fq.append(_get_fastq_files(bc_dir, {'lane':lane.get_name()}, fc_name))
    else:
        for bc in multiplex:
            if not os.path.exists(bc_dir):
                raise IOError("No barcode directory found: " + str(bc_dir))
            item = {'lane':lane.get_name()}
            fq.append(_get_fastq_files(bc_dir, None, item, fc_name, bc_name=bc.get_barcode_id()))
    return fq

def _convert_barcode_id_to_name(multiplex, fc_name, fq):
    bcid2name = dict([(str(mp.get_barcode_id()), mp.get_full_name()) for mp in multiplex])
    if options.no_full_names:
        bcid2name = dict([(str(mp.get_barcode_id()), get_sample_name(mp.get_barcode_name())) for mp in multiplex])
    bcid = re.search("_(\d+)_(\d+)_fastq.txt", fq)
    from_str = "%s_%s_fastq.txt" % (bcid.group(1), bcid.group(2))
    to_str   = "%s_%s.fastq" % (bcid2name[bcid.group(1)], bcid.group(2))
    return fq.replace(from_str, to_str)
 
def _deliver_fastq_file(fq_src, fq_tgt, outdir, fc_link_dir=None):
    f2 = None
    if options.link:
        f = os.symlink
    elif options.move:
        f = shutil.move
    elif options.move_and_relink:
        f = shutil.move
        f2 = os.symlink
    else:
        f = shutil.copyfile
    _handle_data(fq_src, os.path.join(outdir, fq_tgt), f, f2)

def _make_delivery_directory(fc):
    """Make the output directory"""
    _make_dir(fc.get_fc_dir(), "flowcell delivery")
    if (fc.get_fc_alias() != fc.get_fc_id()):
        _handle_data(fc.get_fc_dir(), os.path.join(os.path.dirname(fc.get_fc_dir()), fc.get_fc_alias()), os.symlink)     

def _make_dir(dir, label):
    if not os.path.exists(dir):
        if not options.dry_run:
            os.makedirs(dir)
        logger.info("Creating %s directory %s" % (label, dir))
    else:
        logger.warn("%s already exists: not creating new directory" % (dir))

def _handle_data(src, tgt, f=shutil.copyfile, f2=None):
    if options.only_run_info:
        return
    if src is None:
        return
    if os.path.exists(tgt):
        logger.warn("%s already exists: not doing anything!" %(tgt))
        return
    if options.dry_run:
        print "DRY_RUN: %s file %s to %s" % (f.__name__, src, tgt)
        if not f2 is None:
            print "DRY_RUN: %s file %s to %s" % (f2.__name__, tgt, src)
    else:
        logger.info("%s file %s to %s" % (f.__name__, src, tgt))
        f(src, tgt)
        if not f2 is None:
            logger.info("%s file %s to %s" % (f2.__name__, tgt, src))
            f2(tgt, src)
            

def _deliver_data(data, fastqc, outdir):
    """Loop over data and fastqc and deliver files"""
    for src in data:
        tgt = os.path.join(outdir, os.path.basename(src))
        _handle_data(src, tgt, f=shutil.move if options.move else shutil.copyfile)

    for src in fastqc:
        tgt = os.path.join(outdir, "fastqc", os.path.basename(src))
        _handle_data(src, tgt, f=shutil.move if options.move else shutil.copytree)

def _get_analysis_results(fc, lane):
    """Get analysis results
    
    For now just glob the analysis directory for fastqc output and files with the give flowcell name
    """
    flowcell = "_".join([lane.get_name(), fc.get_fc_id()])
    glob_str = os.path.join(fc.get_fc_dir(), flowcell + "*.*")
    data = glob.glob(glob_str)
    glob_str = os.path.join(fc.get_fc_dir(), "fastqc", flowcell + "*")
    fastqc = glob.glob(glob_str)
    return data, fastqc

def _save_run_info(fc, filename):
    outfile = os.path.join(fc.get_fc_dir(), filename)
    if not options.dry_run:
        with open(outfile, "w") as out_handle:
            yaml.dump(fc.to_structure(), stream=out_handle)
    else:
        print "DRY_RUN:"
        yaml.dump(fc.to_structure(), stream=sys.stdout)

def _sample_based_run_info(fc):
    """Convert a regular run info file into a sample-based info file, one sample per lane"""
    sample_fc = PostProcessedFlowcell(fc.get_fc_name(), fc.get_fc_date(), {}, fc_alias = fc.get_fc_alias(), fc_dir=fc.get_fc_dir())    
    lane_num = 0
    for l in fc.get_lanes():
        bcids = l.get_barcode_ids()
        for barcode_id in bcids:
            s = l.get_sample_by_barcode(barcode_id)
            lane_num = lane_num + 1
            newl = Lane(data={"description":s.get_name(), "lane" :lane_num, "multiplex":[], "analysis":"Minimal", "genome_build":s.get_genome_build(), "sample_prj":s.get_project()})
            newl.set_name("%s" % lane_num)
            files = l.get_files()
            if options.customer_delivery or options.barcode_id_to_name:
                pat = "%s.*_%s_[12].*$" % (l.get_name(), get_sample_name(s.get_barcode_name()))
            else:
                pat = "%s.*_%s_[12].*$" % (l.get_name(), barcode_id)
            bc_files = re.compile(pat, re.IGNORECASE)
            files = filter(bc_files.search, files) 
            newl.set_files(files)
            sample_fc.add_lane(newl)
    _save_run_info(sample_fc, "sample_project_run_info.yaml")
        

def _get_fastq_files(directory, work_dir, item, fc_name, bc_name=None, glob_ext="_fastq.txt",
                    config=None):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if item.has_key("files") and bc_name is None:
        names = item["files"]
        if isinstance(names, basestring):
            names = [names]
        files = [x if os.path.isabs(x) else os.path.join(directory, x) for x in names]
    else:
        assert fc_name is not None
        lane = item["lane"]
        if bc_name:
            glob_str = "%s_*%s_nophix_%s_*%s" % (lane, fc_name, bc_name, glob_ext)
        else:
            glob_str = "%s_*%s_nophix_*%s" % (lane, fc_name, glob_ext)
        files = glob.glob(os.path.join(directory, glob_str))
        files.sort()
        if len(files) > 2 or len(files) == 0:
            raise ValueError("Did not find correct files for %s %s %s %s" %
                    (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            # TODO: Parallelize using pgzip
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        elif fname.endswith(".bam"):
            ready_files = convert_bam_to_fastq(fname, work_dir, config)
        else:
            assert os.path.exists(fname), fname
            ready_files.append(fname)
    ready_files = [x for x in ready_files if x is not None]
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)


if __name__ == "__main__":
    usage = """
    project_management.py   <YAML config file> <flowcell_dir> <project_dir>
                            [<YAML run information>
                             --flowcell_alias=<flowcell_alias>
                             --project_desc=<project_desc> --symlink
                             --move_data --move_and_relink --only_install_run_info
                             --install_data --no_full_names --dry_run --verbose]

    For more extensive help type project_management.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--flowcell_alias", dest="fc_alias")
    parser.add_option("-y", "--project_desc", dest="project_desc")
    parser.add_option("-d", "--install_data", dest="install_data", action="store_true",
                      default=False)
    parser.add_option("-c", "--customer_delivery", dest="customer_delivery", action="store_true",
                      default=False)
    parser.add_option("-b", "--barcode_id_to_name", dest="barcode_id_to_name", action="store_true",
                      default=False)
    parser.add_option("-i", "--only_install_run_info", dest="only_run_info", action="store_true",
                      default=False)
    parser.add_option("-f", "--no_full_names", dest="no_full_names", action="store_true",
                      default=False)
    parser.add_option("-m", "--move_data", dest="move", action="store_true",
                      default=False)
    parser.add_option("-M", "--move_and_relink", dest="move_and_relink", action="store_true",
                      default=False)
    parser.add_option("-l", "--symlink", dest="link", action="store_true",
                      default=False)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print __doc__
        sys.exit()
    
    kwargs = dict(
        fc_alias = options.fc_alias,
        project_desc = options.project_desc,
        )
    main(*args, **kwargs)
