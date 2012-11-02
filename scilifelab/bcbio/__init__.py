"""bcbio init"""
import copy
from scilifelab.utils.misc import opt_to_dict

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
