#!/usr/bin/env python
"""Provide SNP and indel calling.

This script works on delivered project data. In a project 
folder sequence data is stored as fastq files in the data 
directory, grouped by flow cell like names. Optionally, 
it is possible to provide a directory containing bam 
files. 

Usage:
    exome_pipeline.py <exome pipeline YAML config file> 
                      <pruned YAML run information>
                      <project_name> <fc_dir>


The exome pipeline is similar to the post_processing.yaml file in that it holds
information about programs, algorithms and analyses. The contents are more loosely
defined than for post_processing.

The pruned YAML run information file contains run_info about the samples in this 
project only, and also contains the necessary information for conversion between 
barcode ids and barcode names, should the lanes have been multiplexed.

LANE_DATE_FC[_BCI][_1/2]_fastq.txt -> LANE_DATE_FC[_RUNINFONAME][_1/2].fastq

In general, data is delivered as fastq files to a project. Generally, we adopt 
the following convention:

j_doe_00_01/data/fastq_dir

where j_doe_00_01 is the project name to be provided at the command line. Data 
analyses are then performed in 

j_doe_00_01/intermediate/fastq_dir/

which in turn contains subdirectories for flowcells, alignments etc.

The <fc_dir> names a flowcell directory:

j_doe_00_01/intermediate/fastq_dir/fc_dir

in which the delivered fastq files have been renamed to their original names
and link back to the files in j_doe_00_01/data/fastq_dir directory. Hence,
the bcbio modules can be directly applied to the file names in fc_dir.

Relinking of files is done in the script setup_project_files.py.

Requires:
  - bcftools
  - GATK
  - annovar
"""


import os
import sys
from optparse import OptionParser

import yaml

import subprocess
import glob
import collections

from scilife.log import create_log_handler
from scilife.pipeline import log
from scilife.pipeline import sample
from scilife.pipeline.lane import make_lane_items
from scilife.pipeline.fastq import (map_fastq_barcode, get_samples_from_fastq_dir)

from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio import utils
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.broad import BroadRunner
from bcbio.ngsalign import bwa
from bcbio.pipeline import lane

def main(config_file, fastq_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        run_main(config, config_file, fastq_dir, run_info_yaml)

def run_main(config, config_file, fastq_dir, run_info_yaml):
    work_dir = os.getcwd()
    (_, fq_name) = os.path.split(fastq_dir)

    align_dir = os.path.join(work_dir, "intermediate",  fq_name, "alignments")
    (fc_name, fc_date, run_items) = get_samples_from_fastq_dir(fastq_dir)

    run_info = _get_run_info(fc_name, fc_date, config, run_info_yaml)
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(fastq_dir, config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir, "align": align_dir,
            "work": work_dir, "config": config_dir, "flowcell" : None}

    # Since demultiplexing is already done, just extract run_items
    run_items = run_info['details']
    #(dirs, run_items) = map_fastq_barcode(dirs, run_items)
    
    # Align samples
    # lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    # (
    #     fastq1: '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output/3_110106_FC70BUKAAXX_barcode/3_110106_FC70BUKAAXX_trim_1_fastq.txt', 
    #     fastq2: None, 
    #     genome_build: 'hg19', 
    #     mlane_name: '3_110106_FC70BUKAAXX_trim', 
    #     msample : 'Test3', 
    #     dirs: {'fastq': '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output/../data/automated/../110106_FC70BUKAAXX', 'work': '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output', 'flowcell': '../data/automated/../110106_FC70BUKAAXX', 'config': '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output/../data/automated', 'align': '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output/alignments', 'galaxy': '/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/test_automated_output/../data/automated'},
    #     config : {'galaxy_config': 'universe_wsgi.ini', 'algorithm': {'save_diskspace': True, 'recalibrate': False, 'bc_position': 3, 'num_cores': 1, 'aligner': 'bowtie', 'platform': 'illumina', 'max_errors': 2, 'java_memory': '1g', 'bc_read': 1, 'snpcall': False, 'bc_mismatch': 2}, 'custom_algorithms': {'Minimal': {'aligner': ''}, 'SNP calling': {'recalibrate': True, 'snpcall': True, 'aligner': 'bwa', 'dbsnp': 'snps/dbSNP132.vcf'}}, 'distributed': {'rabbitmq_vhost': 'bionextgen'}, 'analysis': {'towig_script': 'bam_to_wiggle.py'}, 'program': {'samtools': 'samtools', 'bowtie': 'bowtie', 'picard': '/usr/share/java/picard', 'barcode': 'barcode_sort_trim.py', 'pdflatex': 'pdflatex', 'ucsc_bigwig': 'wigToBigWig', 'bwa': 'bwa', 'fastqc': 'fastqc', 'gatk': '/usr/share/java/gatk', 'snpEff': '/usr/share/java/snpeff'}})
    # lane_items = make_lane_items(dirs, config, run_items)
    # for item in lane_items:
    #     print str(item)
    # _run_parallel("process_alignment", lane_items, dirs, config)

    # print "Going to align samples"
    # align_items = ((info, fc_lane, fc_name, fc_date, dirs, config) for info in run_items)
    # _run_parallel("align_sample", align_items, dirs, config)
    sys.exit()
        

    # Process samples
    print "Going to process samples"
    #sample_files, sample_fastq, sample_info = sample.organize_samples(dirs, fc_name, fc_date, run_items)
    # for s in sample_files:
    #     print str(s)
        
        
    # samples = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    # for s in samples:
    #     print "Sample : " + str(s)
    
        #_run_parallel("align_sample", samples, dirs, config)

def _get_run_info(fc_name, fc_date, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        return dict(details=run_details, run_id="")
    else:
        log.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        return galaxy_api.run_details(fc_name, fc_date)


def _run_parallel(fn_name, items, dirs, config):
    """Process a supplied function: single, multi-processor or distributed.
    """
    parallel = config["algorithm"]["num_cores"]
    if str(parallel).lower() == "messaging":
        runner = messaging.runner(dirs, config)
        return runner(fn_name, items)
    else:
        out = []
        fn = globals()[fn_name]
        with utils.cpmap(int(parallel)) as cpmap:
            for data in cpmap(fn, items):
                if data:
                    out.extend(data)
        return out

@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def align_sample(*args):
    return sample.align_sample(*args)

def run_main2(config, config_file, fastq_dir, run_info_yaml, bam_dir, force):
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(align_sample, 
                       ((i, fastq_dir, fc_name, fc_date, align_dir, config, config_file)
                        for i in run_items)):
            pass
    sample_files, sample_fastq, sample_info = organize_samples(align_dir,
            fastq_dir, work_dir, fc_name, fc_date, run_items)
    with utils.cpmap(config["algorithm"]["num_cores"]) as cpmap:
        for _ in cpmap(process_sample,
                       ((name, sample_fastq[name], sample_info[name],
                         bam_files, work_dir, config, config_file)
                        for name, bam_files in sample_files)):
            pass

# This function works on the bam files for a sample and is an
# alternative starting point for the script
@utils.map_wrap
def process_sample(sample_name, fastq_files, info, bam_files, work_dir, config, config_file):
    """Finalize processing for a sample"""
    genome_build = config.get("genome_build", None)
    (_, galaxy_dir) = _get_full_paths("", config, config_file)
    (_, sam_ref)    = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                     galaxy_dir)
    fastq1, fastq2 = _combine_fastq_files(fastq_files, work_dir)
    sort_bam = merge_bam_files(bam_files, work_dir, config)
    #bam_to_wig(sort_bam, config, config_file)
    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        dbsnp_file = get_dbsnp_file(config, sam_ref)
        gatk_bam = recalibrate_quality(sort_bam, sam_ref,
                                       dbsnp_file, config_file)
        log.info("Analyzing recalibration %s" % str(sample_name))
        analyze_recalibration(gatk_bam, fastq1, fastq2)
    for snpcaller in config["snpcalling"]:
        print snpcaller["snpcaller"]
        if snpcaller["snpcaller"]["run"]:
            print "going to run %s" % snpcaller["snpcaller"]


def _combine_fastq_files(in_files, work_dir):
    if len(in_files) == 1:
        return in_files[0]
    else:
        cur1, cur2 = in_files[0]
        out1 = os.path.join(work_dir, os.path.basename(cur1))
        out2 = os.path.join(work_dir, os.path.basename(cur2)) if cur2 else None
        if not os.path.exists(out1):
            with open(out1, "a") as out_handle:
                for (cur1, _) in in_files:
                    with open(cur1) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        if out2 and not os.path.exists(out2):
            with open(out2, "a") as out_handle:
                for (_, cur2) in in_files:
                    with open(cur2) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        return out1, out2

def recalibrate_quality(sort_bam_file, sam_ref, dbsnp_file, config_file):
    """Recalibrate alignments with GATK and provide pdf summary.
    """
    bam_file = sort_bam_file.replace("-sort.bam", ".bam")
    cl = ["picard_gatk_recalibrate.py", config_file, sam_ref, bam_file]
    if dbsnp_file:
        cl.append(dbsnp_file)
    subprocess.check_call(cl)
    out_files = glob.glob("%s*-gatkrecal.bam" % os.path.splitext(sort_bam_file)[0])
    assert len(out_files) == 1, out_files
    return out_files[0]

def analyze_recalibration(recal_file, fastq1, fastq2):
    """Provide a pdf report of GATK recalibration of scores.
    """
    cl = ["analyze_quality_recal.py", recal_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)

def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        subprocess.check_call(cl)
    return wig_file

def merge_bam_files(bam_files, work_dir, config):
    """Merge multiple BAM files from a sample into a single BAM for processing.
    """
    out_file = os.path.join(work_dir, os.path.basename(bam_files[0]))
    picard = BroadRunner(config["program"]["picard"],
                         max_memory=config["algorithm"].get("java_memory", ""))
    picard.run_fn("picard_merge", bam_files, out_file)
    for b in bam_files:
        utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
    return out_file

def organize_samples(align_dir, fastq_dir, work_dir, fc_name, fc_date, run_items):
    """Organize BAM output files by sample name."""
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    for sample in run_items:
        sample_name = sample["sample"]
        fname = os.path.join(align_dir, "%s-sort.bam" % sample_name)
        if os.path.exists(fname):
            bams_by_sample[sample_name].append(fname)
            sample_info[sample_name] = sample_name
            fq1, fq2 = get_fastq_files(fastq_dir, sample_name, fc_name)
            fastq_by_sample[sample_name].append(get_fastq_files(fastq_dir, sample_name, fc_name))
    return sorted(bams_by_sample.items()), dict(fastq_by_sample), sample_info
    


# This function works on one sample and is mapped to the list of samples in run_items
@utils.map_wrap
def align_sample2(info, fastq_dir, fc_name, fc_date, align_dir, config, config_file):
    """Do alignment of a sample"""
    sample_name = info.get("sample", "")
    if (config["algorithm"].get("include_short_name", True) and 
        info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    fastq_dir, galaxy_dir = _get_full_paths(fastq_dir, config, config_file)
    genome_build=config.get("genome_build", None)
    align_ref, sam_ref = get_genome_ref(genome_build,
                                        config["algorithm"]["aligner"], galaxy_dir)
    fastq1, fastq2 = get_fastq_files(fastq_dir, info['sample'], fc_name)
    if os.path.exists(fastq1) and config["algorithm"]["aligner"]:
        do_alignment(fastq1, fastq2, align_ref, sam_ref,
                     sample_name, align_dir, config, config_file)

def do_alignment(fastq1, fastq2, align_ref, sam_ref, sample_name, align_dir, config, config_file):
    """Align to provided reference genome, returning an aligned SAM file.
    Does not use lane information"""
    align_fns = {"bwa" : bwa.align}
    # Currently only allow bwa
    aligner_to_use = config["algorithm"]["aligner"]
    utils.safe_makedir(align_dir)
    log.info("Aligning sample %s with aligner %s" % (sample_name, aligner_to_use))
    align_fn = align_fns[aligner_to_use]
    sam_file = align_fn(fastq1, fastq2, align_ref, sample_name, align_dir, config)
    log.info("Converting sample %s to sorted BAM file" % sample_name)
    sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, sample_name, config, config_file)

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, sample_name,
                    config, config_file):
    """Convert SAM file to sorted BAM file.
    """
    lane, pu = _get_lane_info_from_fastq(fastq1)
    cl = ["picard_sam_to_bam.py", "--name=%s" % sample_name,
            "--rg=%s-%s" % (sample_name, lane), "--pu=%s" % pu,
            config_file, sam_file, ref_file, fastq1]
    if fastq2:
        cl.append(fastq2)
    subprocess.check_call(cl)
    utils.save_diskspace(sam_file, "SAM converted to BAM", config)

def _get_lane_info_from_fastq(fastq):
    fh = os.path.basename(fastq)
    lane = fh.split("_")[0]
    pu = "_".join(fh.split("_")[0:3])
    return lane, pu
    
def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir

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

def get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file


if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)
