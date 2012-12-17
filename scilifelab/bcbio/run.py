"""bcbio run module"""
import os
import re
import yaml
import glob

import pandas as pd

from scilifelab.utils.misc import filtered_walk, query_yes_no, opt_to_dict
from scilifelab.utils.dry import dry_write, dry_backup, dry_unlink, dry_rmdir, dry_makedir
from scilifelab.log import minimal_logger
from scilifelab.bcbio import sort_sample_config_fastq, update_sample_config, update_pp_platform_args, merge_sample_config, get_sample_analysis

LOG = minimal_logger(__name__)

# The analysis script for running the pipeline in parallell mode (on one node)  
PARALLELL_ANALYSIS_SCRIPT="automated_initial_analysis.py"
# The analysis script for running the pipeline in distributed mode (across multiple nodes/cores)
DISTRIBUTED_ANALYSIS_SCRIPT="distributed_nextgen_pipeline.py"
# If True, will sanitize the run_info.yaml configuration file when running non-CASAVA analysis
PROCESS_YAML = True
# If True, will assign the distributed master process and workers to a separate RabbitMQ queue for each flowcell 
FC_SPECIFIC_AMPQ = True
## Name of merged sample output directory
MERGED_SAMPLE_OUTPUT_DIR = "TOTAL"

def _sample_status(x):
    """Find the status of a sample.
    
    Look for output files: currently only look for project-summary.csv"""
    if os.path.exists(os.path.join(os.path.dirname(x), "project-summary.csv")):
        return "PASS"
    else:
        return "FAIL"

def _group_samples(flist, include_merged=False):
    """Group samples by sample name and flowcell
    
    This function assumes flist consists of bcbb-config.yaml files. It reads each file
    and extracts sample name and flowcell for subsequent grouping. Exclude MERGED_SAMPLE_OUTPUT_DIR from grouping.

    :param flist: list of bcbb-config.yaml files

    :returns: dictionary of samples grouped by name and flowcell
    """
    sample_d = {}
    for f in flist:
        if not include_merged and os.path.dirname(f).endswith(MERGED_SAMPLE_OUTPUT_DIR):
            continue
        with open(f) as fh:
            conf = yaml.load(fh)
        if conf.get("details", [])[0].get("multiplex", []):
            sample_id = conf.get("details", [])[0].get("multiplex", [])[0].get("name", None)
        else:
            sample_id = conf.get("details", [])[0].get("description")
        if not sample_id:
            LOG.warn("No sample_id found in file {}; skipping".format(f))
            continue
        if conf.get("details", [])[0].get("flowcell_id", None):
            fc_id = conf.get("details", [])[0].get("flowcell_id")
        else:
            fc_id = conf.get("fc_name", None)
        if not fc_id:
            LOG.warn("No flowcell_id found in file {}; skipping".format(f))
            continue
        if not sample_id in sample_d.keys():
            sample_d[sample_id] = {fc_id:f}
        else:
            sample_d[sample_id][fc_id] = f
    return sample_d

## FIXME: make pandas data frame of all samples, with info about lane,
## flowcell, date, barcode_id, path, sample name -> makes searching for files much easier
def sample_table(flist):
    """Make a table from bcbb-config yaml files.

    :param flist: file list of config files

    :returns: data frame
    """
    samples = []
    for f in flist:
        path = os.path.dirname(f)
        with open(f) as fh:
            conf = yaml.load(fh)
        runinfo = conf.get("details") if conf.get("details", None) else conf
        for info in runinfo:
            lane = info.get("lane", None)
            fc_name = info.get("flowcell_id", None)
            fc_date = info.get("fc_date", None)
            if info.get("multiplex", None):
                for mp in info.get("multiplex"):
                    barcode_id = mp.get("barcode_id", None)
                    sample = mp.get("name", None)
                    samples.append([sample, lane, barcode_id, fc_name, fc_date, path])
            else:
                barcode_id = None
                sample = info.get("description", None)
                fc_name = conf.get("fc_name", None)
                fc_date = conf.get("fc_date", None)
                samples.append([sample, lane, barcode_id, fc_name, fc_date, path])
    return pd.DataFrame(samples, columns=["sample", "lane", "barcode_id", "fc_name", "fc_date", "path"])
                                  
def get_vcf_files(flist, vcfext="sort-gatkrecal-realign-variants-combined"):
    """Get dictionary of vcf files.

    :param flist: yaml sample config file list

    """
    vcf_d = {}
    samples = sample_table(flist)
    grouped = samples.groupby("sample")
    for name, group in grouped:
        LOG.debug("Getting vcf file for sample {}".format(name))
        if len(group) > 1:
            pattern = "*_{}-{}.vcf*".format(MERGED_SAMPLE_OUTPUT_DIR, vcfext)
            path = os.path.join(os.path.dirname(group["path"].values[0]), MERGED_SAMPLE_OUTPUT_DIR)
        else:
            pattern = "{}_*_*{}*_{}*-{}.vcf*".format(group["lane"].values[0], group["fc_name"].values[0], str(group["barcode_id"].values[0]).rstrip(".0"), vcfext)
            path = group["path"].values[0]
        if not os.path.exists(path):
            LOG.info("No merge information for sample {}; skipping".format(name))
        files = glob.glob(os.path.join(path, pattern))
        def vcf_filt(f):
            return re.search(".vcf$|.vcf.gz$", f) != None
        files = [x for x in files if vcf_filt(x)]
        if len(files) == 0:
            LOG.warn("No matching file found for pattern {} in {}; skipping".format(pattern, path))
            continue
        if len(files) > 1:
            LOG.warn("More than 1 matching file found for pattern {}; only using first".format(pattern))
        vcf_d[name] = files[0]
    return vcf_d

def setup_merged_samples(flist, sample_group_fn=_group_samples, **kw):
    """Setup analysis that merges multiple sample runs.

    :param flist: list of file names, by default *-bcbb-config.yaml files
    :param sample_group_fn: function that groups files into samples and sample runs. The function takes flist as input.

    :returns: updated flist with config files for merged samples
    """
    new_flist = []
    sample_d = sample_group_fn(flist)
    for k, v in sample_d.iteritems():
        if len(v) > 1:
            f = v[v.keys()[0]]
            out_d = os.path.join(os.path.dirname(os.path.dirname(f)), MERGED_SAMPLE_OUTPUT_DIR)
            LOG.info("Sample {} has {} sample runs; setting up merge analysis in {}".format(k, len(v), out_d))
            dry_makedir(out_d, dry_run=False)
            pp = kw.get("post_process") if kw.get("post_process", None) else f.replace("-bcbb-config.yaml", "-post_process.yaml")
            with open(pp) as fh:
                conf = yaml.load(fh)
            conf = update_pp_platform_args(conf, **{'jobname': "{}_total".format(k), 'workdir': out_d, 'output': "{}_total-bcbb.log".format(k) })
            pp_new = os.path.join(out_d, os.path.basename(pp))
            dry_unlink(pp_new, dry_run=kw.get('dry_run', True))
            dry_write(pp_new, yaml.safe_dump(conf, default_flow_style=False, allow_unicode=True, width=1000), dry_run=kw.get('dry_run', True))
            ## Setup merged bcbb-config file
            bcbb_config = merge_sample_config(v.values(), sample=k)
            bcbb_config_file = os.path.join(out_d, os.path.basename(v.values()[0]))
            if not os.path.exists(bcbb_config_file) or kw.get('new_config', False):
                dry_unlink(bcbb_config_file, dry_run=kw.get('dry_run', True))
                dry_write(bcbb_config_file, yaml.safe_dump(bcbb_config, default_flow_style=False, allow_unicode=True, width=1000), dry_run=kw.get('dry_run', True))
            ##new_flist.extend(v.values())
            new_flist.extend([bcbb_config_file])
    return new_flist


def find_samples(path, sample=None, pattern = "-bcbb-config.yaml$", only_failed=False, **kw):
    """Find bcbb config files in a path.

    :param path: path to search in
    :param sample: a specific sample, or a file consisting of -bcbb-config.yaml files
    :param pattern: pattern to search for

    :returns: list of file names
    """
    flist = []
    if sample:
        if os.path.exists(sample):
            with open(sample) as fh:
                flist = [x.rstrip() for x in fh.readlines()]
        else:
            pattern = "{}{}".format(sample, pattern)
    def bcbb_yaml_filter(f):
        return re.search(pattern, f) != None
    if not flist:
        flist = filtered_walk(path, bcbb_yaml_filter, exclude_dirs=kw.get("exclude_dirs", None), include_dirs=kw.get("include_dirs", None))
    if only_failed:
        status = {x:_sample_status(x) for x in flist}
        flist = [x for x in flist if _sample_status(x)=="FAIL"]
    if len(flist) == 0 and sample:
        LOG.info("No such sample {}".format(sample))
    return [os.path.abspath(f) for f in flist]

def setup_sample(f, analysis, amplicon=False, genome_build="hg19", **kw):
    """Setup config files, making backups and writing new files

    :param path: root path in which to search for samples
    :param dry_run: dry run flag
    """
    if not os.path.exists(f):
        return
    with open(f) as fh:
        config = yaml.load(fh)
    ## Check for correctly formatted config
    if not config.get("details", None):
        LOG.warn("Couldn't find 'details' section in config file: aborting setup!")
        return

    ## Save file to backup if backup doesn't exist
    f_bak = f.replace("-bcbb-config.yaml", "-bcbb-config.yaml.bak")
    if not os.path.exists(f_bak):
        LOG.info("Making backup of {} in {}".format(f, f_bak))
        dry_backup(os.path.abspath(f), dry_run=kw['dry_run'])

    ## Save command file to backup if it doesn't exist
    cmdf = f.replace("-bcbb-config.yaml", "-bcbb-command.txt")
    if os.path.exists(cmdf):
        cmdf_bak = cmdf.replace("-bcbb-command.txt", "-bcbb-command.txt.bak")
        if not os.path.exists(cmdf_bak):
            LOG.info("Making backup of {} in {}".format(cmdf, cmdf_bak))
            dry_backup(os.path.abspath(cmdf), dry_run=kw['dry_run'])

    ## Save post_process file to backup if it doesn't exist
    ppf = f.replace("-bcbb-config.yaml", "-post_process.yaml")
    if os.path.exists(ppf):
        ppf_bak = ppf.replace("-post_process.yaml", "-post_process.yaml.bak")
        if not os.path.exists(ppf_bak):
            LOG.info("Making backup of {} in {}".format(ppf, ppf_bak))
            dry_backup(ppf, dry_run=kw['dry_run'])

    if analysis:
        config = update_sample_config(config, "analysis", analysis)
    if genome_build:
        config = update_sample_config(config, "genome_build", genome_build)
    config = sort_sample_config_fastq(config)

    ## Remove config file and rewrite
    dry_unlink(f, kw['dry_run'])
    dry_write(f, yaml.safe_dump(config, default_flow_style=False, allow_unicode=True, width=1000), dry_run=kw['dry_run'])

    ## Setup post process only if not provided at command line
    if not kw.get("post_process", None):
        ppfile = f.replace("-bcbb-config.yaml", "-post_process.yaml")
        with open(ppfile) as fh:
            pp = yaml.load(fh)
        ## Need to set working directory to path of bcbb-config.yaml file
        if pp.get('distributed', {}).get('platform_args', None):
            platform_args = pp['distributed']['platform_args'].split()
            if "-D" in platform_args:
                platform_args[platform_args.index("-D")+1] = os.path.dirname(f)
            elif "--workdir" in platform_args:
                platform_args[platform_args.index("--workdir")+1] = os.path.dirname(f)
            pp['distributed']['platform_args'] = " ".join(platform_args)
        ## Change keys for all analyses
        for anl in pp.get('custom_algorithms',{}).keys():
            if kw.get('baits', None):
                pp['custom_algorithms'][anl]['hybrid_bait'] = kw['baits']
            if kw.get('targets', None):
                pp['custom_algorithms'][anl]['hybrid_target'] = kw['targets']
            if amplicon:
                pp['custom_algorithms'][anl]['mark_duplicates'] = False
        if amplicon:
            LOG.info("setting amplicon analysis")
            pp['algorithm']['mark_duplicates'] = False
        if kw.get('galaxy_config', None):
            pp['galaxy_config'] = kw['galaxy_config']
        if kw.get('distributed', None):
            LOG.info("setting distributed execution")
            pp['algorithm']['num_cores'] = 'messaging'
        else:
            LOG.info("setting parallell execution")
            pp['algorithm']['num_cores'] = kw['num_cores']
        if kw.get('snpEff', None):
            LOG.info("setting snpEff to {}".format(kw["snpEff"]))
            pp['program']['snpEff'] = kw['snpEff']
        dry_unlink(ppfile, dry_run=kw['dry_run'])
        dry_write(ppfile, yaml.safe_dump(pp, default_flow_style=False, allow_unicode=True, width=1000), dry_run=kw['dry_run'])

def remove_files(f, **kw):
    ## Remove old files if requested
    keep_files = ["-post_process.yaml$", "-post_process.yaml.bak$", "-bcbb-config.yaml$", "-bcbb-config.yaml.bak$",  "-bcbb-command.txt$", "-bcbb-command.txt.bak$", "_[0-9]+.fastq$", "_[0-9]+.fastq.gz$", "_[0-9]+_fastq.txt.gz$", "_[0-9]+_fastq.txt$",
                  "^[0-9][0-9]_.*.txt$", "JOBID", "PID"]
    pattern = "|".join(keep_files)
    def remove_filter_fn(f):
        return re.search(pattern, f) == None

    workdir = os.path.dirname(f)
    remove_files = filtered_walk(workdir, remove_filter_fn)
    remove_dirs = filtered_walk(workdir, remove_filter_fn, get_dirs=True)
    if len(remove_files) == 0:
        pass
    if len(remove_files) > 0 and query_yes_no("Going to remove {} files and {} directories... Are you sure you want to continue?".format(len(remove_files), len(remove_dirs)), force=kw['force']):
        [dry_unlink(x, dry_run=kw['dry_run']) for x in remove_files]
        ## Sort directories by length so we don't accidentally try to remove a non-empty dir
        [dry_rmdir(x, dry_run=kw['dry_run']) for x in sorted(remove_dirs, key=lambda x: len(x), reverse=True)]

def run_bcbb_command(run_info, post_process=None, **kw):
    """Setup bcbb command to run
    
    :param run_info: run info file 
    :param post_process: post process file
    :param kw: keyword arguments
    
    :returns: command line to run
    """
    if not post_process:
        post_process = run_info.replace("-bcbb-config.yaml", "-post_process.yaml")
    with open(post_process, "r") as fh:
        config = yaml.load(fh)
    if str(config["algorithm"]["num_cores"]) == "messaging":
        analysis_script = DISTRIBUTED_ANALYSIS_SCRIPT
    else:
        analysis_script = PARALLELL_ANALYSIS_SCRIPT
    platform_args = config["distributed"]["platform_args"].split()
    cl = [analysis_script, post_process, os.path.dirname(run_info), run_info]
    return (cl, platform_args)
