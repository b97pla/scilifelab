"""Analysis tools for samples"""

import os
import glob
import collections
from scilife.pipeline.fastq import get_fastq_files
from scilife.pipeline import log
from bcbio.pipeline.alignment import align_to_sort_bam

def align_sample(info, fc_lane, fc_name, fc_date, dirs, config):
    """Do alignment of a sample"""
    sample_name = info.get("sample", "")
    if (config["algorithm"].get("include_short_name", True) and 
        info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)

    genome_build=config.get("genome_build", None)

    log.info("Aligning sample: %s; reference genome: %s; " \
                 "researcher %s; analysis method %s" % 
             (sample_name, genome_build, 
              info.get("researcher", ""), info.get("analysis", "")))

    # align_ref, sam_ref = get_genome_ref(genome_build,
    #                                     config["algorithm"]["aligner"], galaxy_dir)
    fastq1, fastq2 = get_fastq_files(dirs["fastq"], info['sample'], fc_name)

    aligner = config["algorithm"].get("aligner", None)
    if os.path.exists(fastq1) and aligner:
        align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                           fc_lane, sample_name, dirs, config)

def organize_samples(dirs, fc_name, fc_date, run_items):
    """Organize BAM output files by sample name."""
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    for sample in run_items:
        print "sample " + str(sample)
        sample_name = sample["sample"]
        fname = os.path.join(dirs["align"], "%s-sort.bam" % sample_name)
        if os.path.exists(fname):
            bams_by_sample[sample_name].append(fname)
            sample_info[sample_name] = sample_name
            fq1, fq2 = get_fastq_files(dirs["fastq"], sample_name, fc_name)
            fastq_by_sample[sample_name].append(get_fastq_files(dirs["fastq"], sample_name, fc_name))
    return sorted(bams_by_sample.items()), dict(fastq_by_sample), sample_info
    



