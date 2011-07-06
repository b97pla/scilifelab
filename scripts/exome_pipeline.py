#!/usr/bin/env python
"""Provide SNP and indel calling.

This script works on delivered project data. In a project 
folder sequence data is stored as fastq files in the data 
directory, grouped by flow cell like names. Optionally, 
it is possible to provide a directory containing bam 
files. 

Usage:
    exome_pipeline.py <config file> <fastq dir> [bam dir]
    
Requires:
  - bcftools
  - GATK
  - annovar
"""


import os
import sys
import subprocess
import glob
import collections
from optparse import OptionParser
import yaml
import logbook

from bcbio.solexa.flowcell import get_flowcell_info
from bcbio import utils
from bcbio.log import create_log_handler
from bcbio.ngsalign import bwa

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(config_file, fastq_dir, run_info_yaml=None, bam_dir=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
        log_handler = create_log_handler(config, LOG_NAME)
    with log_handler.applicationbound():
        log.debug(config)
        run_main(config, config_file, fastq_dir, run_info_yaml, bam_dir)

def run_main(config, config_file, fastq_dir, run_info_yaml, bam_dir):
    work_dir = os.getcwd()
    (_, fq_name) = os.path.split(fastq_dir)
    if run_info_yaml and os.path.exists(run_info_yaml):
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        run_info = dict(details=run_details, run_id="")
    align_dir = os.path.join(work_dir, fq_name, "alignments")
    # Get fc name and date, and make a dictionary reminiscient of run_items
    # in automated_initial_analysis.py
    (fc_name, fc_date, run_items) = get_fastq_samples(fastq_dir)
    with utils.cpmap(1) as cpmap:
        for _ in cpmap(process_sample, 
                       ((i, fastq_dir, fc_name, fc_date, align_dir, config, config_file)
                        for i in run_items)):
            pass

def get_fastq_samples(fastq_dir, info=None):
    run_items = list()
    if info:
        pass
    else:
        fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq*"))
        if not fastq_files:
            print "No fastq files found in " + fastq_dir
            sys.exit()
        (name, date) = get_flowcell_info(fastq_files[0])
    for fq in fastq_files:
        (_, fc) = os.path.split(fq)
        # files have to be lane_date_flowcell_sample_1/2.fastq
        sample = fc.split("_")[3]
        run_items.append({'sample':sample})

    return name, date, run_items

@utils.map_wrap
def process_sample(info, fastq_dir, fc_name, fc_date, align_dir, config, config_file):
    """Do alignment of a sample"""
    sample_name = info.get("sample", "")
    if (config["algorithm"].get("include_short_name", True) and 
        info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    fastq_dir, galaxy_dir = _get_full_paths(fastq_dir, config, config_file)
    genome_build="hg19"
    align_ref, sam_ref = get_genome_ref(genome_build,
                                        config["algorithm"]["aligner"], galaxy_dir)
    fq1, fq2 = get_fastq_files(fastq_dir, 4, fc_name)
    # print fq1
    # print fq2
    # Process samples
    # In this setting get fastq files and loop over them 
    # bam_sample_files, sample_fastq, sample_info = organize_samples(align_dir,
    #                                                                fastq_dir, work_dir, fc_name, fc_date, run_items)
    

def _add_full_path(dirname, basedir=None):
    if basedir is None:
        basedir = os.getcwd()
    if not dirname.startswith("/"):
        dirname = os.path.join(basedir, dirname)
    return dirname

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = _add_full_path(fastq_dir)
    config_dir = _add_full_path(os.path.dirname(config_file))
    galaxy_config_file = _add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file)

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    if not aligner or not genome_build:
        return (None, None)
    ref_files = dict(
            bowtie = "bowtie_indices.loc",
            bwa = "bwa_index.loc",
            samtools = "sam_fa_indices.loc",
            maq = "bowtie_indices.loc",
            tophat = "bowtie_indices.loc")
    remap_fns = dict(
        maq = _remap_to_maq
        )
    out_info = []
    ref_dir = os.path.join(galaxy_base, "tool-data")
    for ref_get in [aligner, "samtools"]:
        ref_file = os.path.join(ref_dir, ref_files[ref_get])
        cur_ref = None
        with open(ref_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        cur_ref = parts[-1]
                        break
        if cur_ref is None:
            raise IndexError("Genome %s not found in %s" % (genome_build,
                ref_file))
        try:
            cur_ref = remap_fns[ref_get](cur_ref)
        except KeyError:
            pass
        out_info.append(_add_full_path(cur_ref, ref_dir))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)


def _remap_to_maq(ref_file):
    """ToDo: Why is this needed for ?
    """
    base_dir = os.path.dirname(os.path.dirname(ref_file))
    name = os.path.basename(ref_file)
    for ext in ["fa", "fasta"]:
        test_file = os.path.join(base_dir, "maq", "%s.%s" % (name, ext))
        if os.path.exists(test_file):
            return test_file
    raise ValueError("Did not find maq file %s" % ref_file)

def get_sample_info():
    """Get the sample info. For now looks in fastq_dir and groups fastq files according to sample name"""
    pass


def organize_samples(align_dir, fastq_dir, work_dir, fc_name, fc_date, run_items):
    """Organize BAM output files by sample name"""
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    glob_str = "*.fastq"
    #for fq_file in run_items:
        
    


def get_fastq_file2s(fastq_dir, lane, fc_name, bc_name=None):
    pass


if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)
