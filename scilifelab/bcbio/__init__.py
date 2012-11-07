"""bcbio init"""
import os
import copy
from scilifelab.utils.misc import opt_to_dict
from scilifelab.log import minimal_logger

LOG = minimal_logger(__name__)


POST_PROCESS_OPTS = ["-t", "--time", "-A", "--account", "-p", "--partition", "-D", "--workdir", "-o", "--output", "-J", "--job-name"]

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
        for j in range(0, len(runinfo[i].get("multiplex"))):
            LOG.info("Setting {} to {} for sample {}".format(key, value, runinfo[i]["multiplex"][j]["name"]))
            runinfo[i]["multiplex"][j][key] = value
    return {'details':runinfo}
    
def sort_sample_config_fastq(conf):
    """Sort fastq entires in sample config file, at the same
    time adding/subtracting gz extension if present/absent.

    :param conf: sample configuration
    
    :returns: updated post process file
    """
    newconf = copy.deepcopy(conf)
    runinfo = newconf.get("details") if newconf.get("details", None) else newconf
    for i in range(0, len(runinfo)):
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
    return {'details':runinfo}
