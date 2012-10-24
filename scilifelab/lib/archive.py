"""scilifelab lib module"""

import os
import re
import subprocess

from cStringIO import StringIO
from scilifelab.utils.misc import filtered_walk

import scilifelab.log

LOG = scilifelab.log.minimal_logger(__name__)


def flowcell_remove_status(archive_dir, swestore_dir, to_remove="to_remove"):
    """This function looks for flowcells that could be deleted
    from archive and returns a list of flowcells with a KEEP/RM
    flag. The rules are
    
    1. the flowcell is in archive to_remove file
    2. pbzip ran without error
    3. the tarball filesize looks ok
    4. checksum irods is ok

    :param archive_dir: archive directory
    :param swestore_dir: base dir for swestore
    :param to_remove: to remove file name
    """
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}
    ## Check for ils
    try:
        proc = subprocess.Popen(["ils"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        proc.wait()
        proc = subprocess.Popen(["icd", os.path.basename(os.path.dirname(archive_dir))],  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        proc.wait()
    except:
        LOG.warn("No such command 'ils': please load the irods module" )
        return output_data
    ## make flowcell dictionary based on to_remove contents
    to_remove_file = os.path.join(archive_dir, to_remove)
    with open(to_remove_file) as fh:
        remove_list = fh.readlines()
    flowcells = {k.replace("./", "").rstrip():{'in_archive':False, 'pbzip_exit':1, 'tarball_size':0, 'irods_checksum':1} for k in remove_list if k.rstrip() != ''}

    ## Look for compress logs
    pattern = "slurm.*.out$"
    def compress_fn(f):
        return re.search(pattern, f) != None
    compress_log_files = filtered_walk(os.path.join(archive_dir, "compress_logs"), compress_fn)
    for f in compress_log_files:
        with open(f) as fh:
            compress_str = "".join([x.strip() for x in fh.readlines()])
        m = re.search("Compressing[ ]+([0-9A-Za-z_\-]+)\.\.\..*Exit code:[ ]+([0-9]+)", compress_str)
        if m:
            if not m.groups()[0] in flowcells.keys():
                LOG.warn("flowcell {} present in to_remove but not in archive".format(m.groups()[0]))
            else:
                flowcells[m.groups()[0]]['pbzip_exit'] = m.groups()[1]
        else:
            LOG.warn("{}: no match for {}".format(f, compress_str))

    ## Get tarball sizes and check if in archive 
    ## Loop through flowcells and perform ichksum
    for k in flowcells.keys():
        LOG.debug("Getting tarball size, archive presence and ichksum for {}".format(k))
        fcdir = os.path.join(archive_dir, k)
        if os.path.exists(fcdir):
            flowcells[k]['in_archive'] = True
        fctar = os.path.join(swestore_dir, "drophere2archive", "{}.tar.bz2".format(k))
        try:
            cl = ["ichksum", os.path.basename(fctar)]
            proc = subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, stderr) = proc.communicate()
            proc.wait()
            flowcells[k]['irods_checksum'] = stdout.split("\n")[1]
        except:
            LOG.warn("command {} failed".format(" ".join(cl)))
        if not os.path.exists(fctar):
            continue
        else:
            LOG.debug("tarball exists: {}".format(fctar))
            statinfo = os.stat(fctar)
            flowcells[k]['tarball_size'] = float(int(statinfo.st_size) / 1e9)

    output_data["stdout"].write("\nFlowcell archive status\n")
    output_data["stdout"].write("=======================\n")
    output_data["stdout"].write("\nThe table lists those flowcells still present in archive. The exict code for pbzip should be 0\nfor success. A non-existing tarball has size 0.\n\n")
    output_data["stdout"].write("{:<40}{:>12}{:>20}{:>60}\n".format("Flowcell", "pbzip_exit", "tarball_size (G)", 'irods_checksum'))
    output_data["stdout"].write("{:<40}{:>12}{:>20}{:>60}\n".format("--------", "----------", "----------------", '--------------'))

    for k in sorted(flowcells.keys()):
        if not flowcells[k]['in_archive']:
            continue
        output_data["stdout"].write("{:<40}{:>12}{:>20.2f}{:>60}\n".format(k, flowcells[k]['pbzip_exit'], flowcells[k]['tarball_size'], flowcells[k]['irods_checksum'] ))
    return output_data

