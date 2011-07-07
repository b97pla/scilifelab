#!/usr/bin/env python
"""Provide SNP and indel calling.

This script works on delivered project data. In a project 
folder sequence data is stored as fastq files in the data 
directory, grouped by flow cell like names. Optionally, 
it is possible to provide a directory containing bam 
files. 

Usage:
    exome_pipeline.py <config file> <fastq dir> [<bam dir>] 
    [--force]
    
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
from bcbio.broad import BroadRunner
from bcbio.log import create_log_handler
from bcbio.ngsalign import bwa

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(config_file, fastq_dir, run_info_yaml=None, bam_dir=None, force=False):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
        log_handler = create_log_handler(config, LOG_NAME)
    with log_handler.applicationbound():
        #log.debug(config)
        run_main(config, config_file, fastq_dir, run_info_yaml, bam_dir, force)

def run_main(config, config_file, fastq_dir, run_info_yaml, bam_dir, force):
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
    for item in run_items:
        log.info("Going to process sample " + item['sample']) 
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
def align_sample(info, fastq_dir, fc_name, fc_date, align_dir, config, config_file):
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
    
def get_fastq_files(directory, sample, fc_name, bc_name=None):
    """Retrieve fastq files for the given sample, ready to process.
    """
    if bc_name:
        glob_str = "*_*_%s_%s_%s_*.fastq*" % (fc_name, sample, bc_name)
    else:
        glob_str = "*_*_%s_%s_*.fastq*" % (fc_name, sample)
    files = glob.glob(os.path.join(directory, glob_str))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, sample, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

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

def get_fastq_samples(fastq_dir, info=None):
    run_items = list()
    samples = {}
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
        if not samples.has_key(sample):
            run_items.append({'sample':sample})
        samples[sample] = True
    return name, date, run_items

def get_dbsnp_file(config, sam_ref):
    snp_file = config["algorithm"].get("dbsnp", None)
    if snp_file:
        base_dir = os.path.dirname(os.path.dirname(sam_ref))
        snp_file = os.path.join(base_dir, snp_file)
    return snp_file


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--force", dest="force", action="store_true", default=False)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict(
        force=options.force)
    main(*args, **kwargs)
