#!/usr/bin/env python
"""Deliver results from a project based analysis.

Usage:
  project_data_delivery.py <YAML config file> <delivery_dir> <YAML run info>
                            [--analysis_dir --move_data --dry_run --verbose --no_vcf
                             --no_metrics --no_bigwig --no_rename --bam --bam_glob]

Given a YAML config file, a delivery directory, and a YAML run info
file describing analysis setup, copy files to delivery directory. By
default, the following files are delivered:

  - project-summary.csv
  - YAML run info if present
  - run_summary.yaml
  - *.vcf
  - *.vcf.idx
  - *metrics
  - *.bigwig

In addition, bam files can be delivered by using the --bam flag. The
setup of project based analyses link lane numbers to samples, and by
default, the script will rename the filenames from
LANE_ANALYSISDATE_ANALYSISID to LANE_SAMPLEID_ANALYSISDATE_ANALYSISID.
The SAMPLEID is found in the description field of the YAML run info.

Options:
  -a, --analysis_dir=<analysis dir>             Analysis directory where results are found. Default
                                                current working directory.
  -m, --move_data                               Move data instead of copying
  -V, --no_vcf                                  Don't deliver vcf files
  -M, --no_metrics                              Don't deliver metrics files
  -B, --no_bigwig                               Don't deliver bigwig files
  -R, --no_rename                               Don't rename files
  -b, --bam                                     Deliver bam files
  -g, --bam_glob                                Glob to use for bam files.
                                                Default: *-sort-dup-gatkrecal-realign*.bam
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

from bcbio.log import logger, setup_logging
from bcbio.pipeline.config_loader import load_config

def main(config_file, delivery_dir, run_info_yaml, analysis_dir=None):
    if analysis_dir is None:
        analysis_dir = os.path.abspath(os.path.curdir)
    config = load_config(config_file)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(analysis_dir, "log")
    setup_logging(config)

    if not os.path.exists(analysis_dir):
        logger.error("No analysis directory found!")
        sys.exit()

    with open(run_info_yaml) as fp:
        run_info_structure = yaml.load(fp)    
    
    lane2sample = dict()
    infiles = dict()
    for info in run_info_structure['details']:
        lane2sample[info['lane']] = info['description']
        infiles[info['lane']] = dict(vcf=[], bam=[], bigwig=[], metrics=[])
        # Vcf files, tsv and tranches
        vcftypes = ('*.vcf', '*.idx', '*.tranches', '*.eval', '*.tsv')
        for vcftype in vcftypes:
            glob_str = os.path.join(analysis_dir, str(info['lane']) +  "_" + vcftype)
            infiles[info['lane']]['vcf'].extend(glob.glob(glob_str))
        # Bam files
        glob_str = os.path.join(analysis_dir, str(info['lane']) + "_" + options.bam_glob)
        bamfiles = glob.glob(glob_str)
        infiles[info['lane']]['bam'] = bamfiles
        # Bigwig files
        glob_str = os.path.join(analysis_dir, str(info['lane']) + "_" + "*.bigwig")
        bigwigfiles = glob.glob(glob_str)
        infiles[info['lane']]['bigwig'] = bigwigfiles
        # metrics files
        glob_str = os.path.join(analysis_dir, str(info['lane']) +  "_" + "*metrics")
        metricsfiles = glob.glob(glob_str)
        infiles[info['lane']]['metrics'] = metricsfiles

    # snpEff files
    glob_str = os.path.join(analysis_dir, "snpEff*")
    snpeff_files = glob.glob(glob_str)
    
    # Loop through the list and deliver if appropriate
    _make_dir(delivery_dir)
    _deliver_file(os.path.join(analysis_dir,"project-summary.csv"),os.path.join(delivery_dir,"project-summary.csv") )
    _deliver_file(os.path.join(analysis_dir,"run_summary.yaml"),os.path.join(delivery_dir,"run_summary.yaml") )
    _deliver_file(run_info_yaml,os.path.join(delivery_dir,os.path.basename(run_info_yaml) ))
    if not options.no_vcf:
        for sf in snpeff_files:
            _deliver_file(sf, os.path.join(delivery_dir,os.path.basename(sf) ))
    for lane_num in infiles.keys():
        lane = infiles[lane_num]
        if not options.no_vcf:
            for vcf in lane['vcf']:
                (src, tgt) = _rename_sample_file(vcf, lane_num, lane2sample[lane_num],delivery_dir)
                _deliver_file(src, tgt)
        if not options.no_bigwig:
            for bigwig in lane['bigwig']:
                (src, tgt) = _rename_sample_file(bigwig, lane_num, lane2sample[lane_num],delivery_dir)
                _deliver_file(src, tgt)
        if not options.no_metrics:
            for metrics in lane['metrics']:
                (src, tgt) = _rename_sample_file(metrics, lane_num, lane2sample[lane_num],delivery_dir)
                _deliver_file(src, tgt)
        if options.bam:
            for bamfile in lane['bam']:
                (src, tgt) = _rename_sample_file(bamfile, lane_num, lane2sample[lane_num],delivery_dir)
                _deliver_file(src, tgt)



def _rename_sample_file(src, lane, sample, outdir):
    if options.no_rename:
        return (os.path.abspath(src), os.path.abspath(os.path.join(outdir, os.path.basename(src))))
    from_str = "%s_" % (lane)
    to_str = "%s_%s_" % (lane, sample)
    tgt = os.path.basename(src).replace(from_str, to_str, 1)
    return (os.path.abspath(src), os.path.abspath(os.path.join(outdir, tgt)))

def _deliver_file(src, tgt):
    if options.move:
        deliver_fn = shutil.move
    else:
        deliver_fn = shutil.copyfile
    if src is None:
        return
    if not os.path.exists(src):
        return
    if os.path.exists(tgt):
        logger.warn("%s already exists: not doing anything!" %(tgt))
        return
    if options.dry_run:
        print "DRY_RUN: %s file %s to %s" % (deliver_fn.__name__, src, tgt)
    else:
        logger.info("%s file %s to %s" % (deliver_fn.__name__, src, tgt))
        f(src, tgt)


def _make_dir(dir):
    if options.dry_run:
        return
    if not os.path.exists(dir):
        os.makedirs(dir)
        logger.info("Creating delivery directory %s" % (dir))
    else:
        logger.warn("%s already exists: not creating new directory" % (dir))

if __name__ == "__main__":
    usage = """
    project_data_delivery.py <YAML config file> <delivery_dir>
                             [--analysis_dir --move_data --dry_run --verbose --no_vcf
                              --no_metrics --no_bigwig --no_rename --bam --bam_glob]

    For more extensive help type project_data_delivery.py
"""
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--analysis_dir", dest="analysis_dir")
    parser.add_option("-b", "--bam", dest="bam", action="store_true",
                      default=False)
    parser.add_option("-g", "--bam_glob", dest="bam_glob", default="*-sort-dup-gatkrecal-realign*.bam*")
    parser.add_option("-f", "--only_install_run_info", dest="only_run_info", action="store_true",
                      default=False)
    parser.add_option("-m", "--move_data", dest="move", action="store_true",
                      default=False)
    parser.add_option("-V", "--no_vcf", dest="no_vcf", action="store_true",
                      default=False)
    parser.add_option("-M", "--no_metrics", dest="no_metrics", action="store_true",
                      default=False)
    parser.add_option("-B", "--no_bigwig", dest="no_bigwig", action="store_true",
                      default=False)
    parser.add_option("-R", "--no_rename", dest="no_rename", action="store_true",
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
        analysis_dir = options.analysis_dir,
        )
    main(*args, **kwargs)
