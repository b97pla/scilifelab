"""bcbio init"""
import os
import copy
import yaml
from itertools import izip
from scilifelab.utils.dry import dry_rsync
from scilifelab.utils.misc import opt_to_dict
from scilifelab.log import minimal_logger

from bcbio.pipeline.run_info import _unique_flowcell_info

LOG = minimal_logger(__name__)

POST_PROCESS_OPTS = ["-t", "--time", "-A", "--account", "-p", "--partition", "-D", "--workdir", "-o", "--output", "-J", "--job-name"]
POST_PROCESS_OPT_GROUPS = {"time":["-t", "--time"], 
                           "account":["-A", "--account"],
                           "partition":["-p", "--partition"],
                           "workdir":["-D", "--workdir"],
                           "jobname":["-J", "--job-name"],
                           "output":["-o", "--output"]
                           }

def update_pp_platform_args(conf, **kw):
    """Update post process platform args with keys given in kw.

    :param conf: post process configuration
    :param kw: keywords with platform arg keys
    """
    newconf = copy.deepcopy(conf)
    opt_d = opt_to_dict(conf.get("distributed", {}).get("platform_args", []).split())
    for k in kw.keys():
        if not k in POST_PROCESS_OPT_GROUPS.keys():
            LOG.warn("Trying to update disallowed key {}; skipping".format(k))
        else:
            for x in POST_PROCESS_OPT_GROUPS[k]:
                if x in opt_d.keys(): 
                    del opt_d[x]
            opt_d[POST_PROCESS_OPT_GROUPS[k][0]] = kw[k]
    platform_args = ["{}={}".format(x, opt_d[x]) if x.startswith("--") else "{} {}".format(x, opt_d[x]) for x in opt_d.keys()]
    newconf["distributed"]["platform_args"] = " ".join(platform_args)
    return newconf
            
def prune_pp_platform_args(conf, keep_opts = POST_PROCESS_OPTS):
    """Rewrite platform args of a post_process file.

    :param conf: post process config object
    :param keep_opts: options to keep in platform args

    :returns: updated post process object
    """
    newconf = copy.deepcopy(conf)
    if conf.get("distributed", {}).get("platform_args"):
        opt_d = opt_to_dict(conf['distributed']['platform_args'].split())
        keys = list(set(keep_opts) & set(opt_d.keys()))
        platform_args = ["{}={}".format(x, opt_d[x]) if x.startswith("--") else "{} {}".format(x, opt_d[x]) for x in keys]
        newconf['distributed']['platform_args'] = " ".join(platform_args)
    return newconf

def update_sample_config(conf, key, value):
    """Update sample config dictionary fields 'key' to value.

    :param conf: sample configuration. 
    :param key: key to update
    :param value: update key to value

    :returns: updated configuration
    """
    newconf = copy.deepcopy(conf)
    runinfo = newconf.get("details") if newconf.get("details", None) else newconf
    for i in range(0, len(runinfo)):
        runinfo[i][key] = value
        if runinfo[i].get("multiplex", None):
            for j in range(0, len(runinfo[i].get("multiplex"))):
                LOG.info("Setting {} to {} for sample {}".format(key, value, runinfo[i]["multiplex"][j]["name"]))
                runinfo[i]["multiplex"][j][key] = value
    newconf['details'] = runinfo
    return newconf

def get_sample_analysis(conf):
    """Get the analysis defined in a config file.

    :param conf: bcbb config file

    :returns: list of analyses
    """
    analysis = []
    runinfo = conf.get("details") if conf.get("details", None) else conf
    for i in range(0, len(runinfo)):
        analysis.extend([runinfo[i].get("analysis", None)])
        for j in range(0, len(runinfo[i].get("multiplex"))):
            analysis.extend([runinfo[i]["multiplex"][j].get("analysis", None)])
        analysis.extend([None])
    return list(set(analysis) - set([None]))
    
    
def sort_sample_config_fastq(conf):
    """Sort fastq entires in sample config file, at the same
    time adding/subtracting gz extension if present/absent.

    :param conf: sample configuration
    
    :returns: updated post process file
    """
    newconf = copy.deepcopy(conf)
    runinfo = newconf.get("details") if newconf.get("details", None) else newconf
    for i in range(0, len(runinfo)):
        if runinfo[i].get("multiplex", None):
            for j in range(0, len(runinfo[i].get("multiplex"))):
                runinfo[i]["multiplex"][j]["files"].sort()
                seqfiles = runinfo[i]["multiplex"][0]["files"]
                for k in range(0, len(seqfiles)):
                    if not os.path.exists(os.path.abspath(seqfiles[k])):
                        (_, ext) = os.path.splitext(seqfiles[k])
                        if ext == ".gz":
                            LOG.warn("Couldn't find gz file; will set use extension .fastq")
                            runinfo[i]["multiplex"][j]["files"][k] = runinfo[i]["multiplex"][j]["files"][k].replace(".gz", "")
                        else:
                            LOG.warn("Couldn't find fastq file; will set use extension .gz")
                            runinfo[i]["multiplex"][j]["files"][k] = "{}.gz".format(runinfo[i]["multiplex"][j]["files"][k])
        else:
            ## Assume files directly in lane
            seqfiles = runinfo[i]["files"]
            for k in range(0, len(seqfiles)):
                if not os.path.exists(os.path.abspath(seqfiles[k])):
                    (_, ext) = os.path.splitext(seqfiles[k])
                    if ext == ".gz":
                        LOG.warn("Couldn't find gz file; will set use extension .fastq")
                        runinfo[i]["files"][k] = runinfo[i]["files"][k].replace(".gz", "")
                    else:
                        LOG.warn("Couldn't find fastq file; will set use extension .gz")
                        runinfo[i]["files"][k] = "{}.gz".format(runinfo[i]["files"][k])
    newconf['details'] = runinfo
    return newconf

def merge_sample_config(flist, sample, out_d, dry_run=True):
    """Merge sample config files, making unique lanes if necessary.

    Also copies sequence files with rsync to the output directory.
    This is a workaround for the case where sequence file names are
    identical for different flowcell runs, causing the pipeline to
    crash.

    :param flist: list of configuration files
    :param sample: sample name to be used in description field for merging
    :param out_d: output directory
    
    :returns: merged configuration 
    """
    newconf = {'details':[]}
    lane = 1
    for f in flist:
        with open(f) as fh:
            conf = yaml.load(fh)
        # Make sure the fastq files exist
        conf = sort_sample_config_fastq(conf)
        runinfo = conf.get("details") if conf.get("details", None) else conf
        for i in range(0, len(runinfo)):
            for j in range(0, len(runinfo[i].get("multiplex"))):
                seqfiles = [os.path.join(os.path.dirname(f), x) for x in runinfo[i]["multiplex"][0]["files"]]
                target_seqfiles = [os.path.join(out_d, os.path.basename(x).replace(sample, "{}_{}".format(sample, runinfo[i]["flowcell_id"]))) for x in seqfiles]
                [dry_rsync(src, tgt, dry_run=dry_run) for src, tgt in izip(seqfiles, target_seqfiles)]
                info = {}
                info["lane"] = str(lane)
                info["analysis"] = runinfo[i]["analysis"]
                info["description"] = str(sample)
                info["files"] = target_seqfiles
                info["genome_build"] = runinfo[i]["genome_build"]
                newconf['details'].append(info)
                lane = lane + 1
    (fc_name, fc_date) = _unique_flowcell_info()
    newconf['fc_date'] = fc_date
    newconf['fc_name'] = "TOTAL"
    return newconf
