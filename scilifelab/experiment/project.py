"""scilifelab.experiment.project module"""
import os
import re
from scilifelab.utils.misc import filtered_walk
# By default samples are named P###.*
re_sample = "^(P[0-9]{3,})"

def validate_sample_directories(plist, pdir):
    """Validate that paths in plist are indeed in subdirectories of pdir.

    :param plist: list of paths
    :param pdir: project/production directory

    :raises: Exception if item in plist not in subdirectory of pdir
    """
    for path in plist:
        if not os.path.realpath(path).startswith(os.path.realpath(pdir)) or os.path.realpath(path) == os.path.realpath(pdir):
            LOG.warning("Path {} *must* be in a subdirectory of {}".format(path, pdir))
            raise Exception

def find_samples(path, sample=None, pattern=re_sample, depth=1, **kw):
    """Look for samples stored in directories given a path. The
    directory name is assumed to be equal to the sample name.

    :param path: path to search in
    :param sample: a specific sample or a file consisting of sample names
    :param depth: depth to look for samples. NB: currently not implemented
    :param re_sample: regexp to use when looking for samples

    :returns: a list of sample directory paths
    """
    def sample_name_filter(f):
        return re.search(pattern, f) != None
    plist = []
    if sample:
        if os.path.exists(sample):
            with open(sample) as fh:
                samplelist = fh.readlines()
            plist = [os.path.join(path, x.rstrip()) for x in samplelist if re.search(pattern, x)]
            if len(plist) == 0:
                return plist
        else:
            if os.path.exists(os.path.join(path, sample)):
                plist = [os.path.join(path, sample)]
            else:
                return plist
    if not plist:
        potential_samples = os.listdir(path)
        plist = [os.path.join(path, os.path.basename(x)) for x in potential_samples if sample_name_filter(x)]
    validate_sample_directories(plist, path)
    return plist
        
