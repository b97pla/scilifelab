"""Pm module clean"""
import os
import re

from scilifelab.utils.misc import query_yes_no, filtered_walk
from scilifelab.utils.dry import dry_unlink, dry_write

from cement.core import backend
LOG = backend.minimal_logger(__name__)

MINFILESIZE = 2000

def _purge_by_sample(files, dry_run, fsize=MINFILESIZE):
    saved_size = 0
    for i in range(0, len(files)-1):
        f1 = os.path.basename(files[i])
        f2 = os.path.basename(files[i])
        if f1.startswith(os.path.splitext(f2)[0]):
            statinfo = os.stat(files[i])
            if statinfo.st_size < fsize:
                continue
            saved_size = saved_size + statinfo.st_size
            LOG.info("Purging bam file {}".format(files[i])) 
            dry_unlink(files[i], dry_run)
            dry_write(files[i], "File removed to save disk space: Moved to {}".format(files[i+1]), dry_run)
    return saved_size

def purge_alignments(path, ftype="sam", keep="last", dry_run=False, force=False, fsize=MINFILESIZE):
    """Cleanup sam and bam files. In some cases, sam files persist. If
    the corresponding bam file exists, replace the sam file contents
    with a message that the file has been removed to save space.
    
    In general, several bam files are produced in an analysis. By
    grouping bam files by prefix, either the most recent file is
    retained for further reference, or a specific analysis is kept.
    """
    if ftype == "sam":
        pattern = ".sam$"
    elif ftype == "bam":
        pattern = ".bam$"
    else:
        LOG.warn("ftype must be one of 'sam' or 'bam'")
        return
    LOG.debug("running purge_alignments in path {} with pattern {} keep rule {}".format(path, pattern, keep))
    def purge_filter(f):
        if not pattern:
            return
        return re.search(pattern, f) != None
    
    flist = filtered_walk(path, purge_filter, exclude_dirs=["realign-split"])
    if len(flist) == 0:
        LOG.info("No {} files found in {}".format(ftype, path))
        return
    if len(flist) > 0 and not query_yes_no("Going to remove/cleanup {} {} files ({}...). Are you sure you want to continue?".format(len(flist), ftype, ",".join([os.path.basename(x) for x in flist[0:10]])), force=force):
        return
    if ftype == "sam":
        for f in flist:
            LOG.info("Purging {} file {}".format(ftype, f))
            dry_unlink(f, dry_run)
            if os.path.exists(f.replace(".sam", ".bam")):
                dry_write(f, "File removed to save disk space: SAM converted to BAM", dry_run)
        return
    elif ftype == "bam":
        samples = {}
        for f in flist:
            m = re.search("([0-9A-Za-z\_]+)-.*", os.path.basename(f))
            if not m:
                LOG.debug("Couldn't determine prefix for {}".format(f))
                continue
            sid = m.groups()[0]
            
            if not sid in samples.keys():
                samples[sid] = {}
            dname = os.path.dirname(f) 
            if not dname in samples[sid].keys():
                samples[sid][dname] = []
            samples[sid][dname].append(f)

        saved_size = 0
        for k in samples.iterkeys():
            for d, files  in samples[k].iteritems():
                if not files or len(files) == 1:
                    continue
                files.sort(lambda x,y: cmp(len(x), len(y)))
                if keep == "last":
                    LOG.info("Keeping file {} and removing all files with common prefix: {}".format(os.path.basename(files[len(files)-1]), ", ".join([os.path.basename(x) for x in files[0:-1]])))
                saved_size = _purge_by_sample(files, dry_run, int(fsize)) + saved_size
        LOG.info("Will save approximately {:.1f}G space".format(saved_size / 1e9))
