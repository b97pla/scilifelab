"""scilifelab lib module"""

import os
import re
import subprocess
import shlex

from cStringIO import StringIO
from fabric.api import task, run, execute, cd, settings, local
from fabric.network import disconnect_all
from datetime import datetime

from scilifelab.utils.misc import filtered_walk
from scilifelab.utils.misc import query_yes_no, md5sum

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

def rm_run(arch, root, flowcell=None):
    """Remove a flowcell folder from the root folder
    """
    path = os.path.join(root,flowcell)
    if not query_yes_no("Going to remove flowcell folder {}. This action can not be undone. Are you sure you want to continue?".format(path),
                        force=arch.pargs.force):
        return
    arch.log.info("removing {}".format(path))
    arch.app.cmd.rmtree(path)

def rm_tarball(arch, tarball):
    """Remove a tarball
    """
    if not query_yes_no("Going to remove tarball {}. This action can not be undone. Are you sure you want to continue?".format(tarball),
                        force=arch.pargs.force):
        return
    arch.log.info("removing {}".format(tarball))
    arch.app.cmd.safe_unlink(tarball)

def package_run(arch, root, flowcell, workdir=None, excludes=None, compress_program=None, check_finished=False, **kw):
    """Package a run in preparation for archiving to swestore
    """

    # Check that a supplied file with excludes exist
    if excludes is not None and not os.path.exists(excludes):
        arch.log.error("Excludes file {} does not exist".format(excludes))
        return None

    if not workdir:
        workdir = root
    elif not os.path.exists(workdir):
        arch.log.info("Creating non-existing working directory {}".format(workdir))
        arch.app.cmd.safe_makedir(workdir)

    if check_finished:
        if not os.path.exists(os.path.join(root, flowcell, 'RTAComplete.txt')):
            arch.log.error("check_finished option was enabled and the run has not " \
                "finished synching yet.")
            return None
        else:
            #Even if RTAComplete.txt is present we give a threashold of 24 hours
            #just in case we're transfering the whole run. We cannot ensure that
            #RTAComplete.txt is the last file to be transfered
            stats = os.stat(os.path.join(root, flowcell, 'RTAComplete.txt'))
            mod_time = datetime.now() - datetime.fromtimestamp(stats.st_mtime)
            if not mod_time.days:
                arch.log.warn("RTAComplete.txt file is present, but it is not " \
                              "older than 1 day. Not packaging the run.")
                return None

    dest_path = os.path.join(workdir,"{}.tar{}".format(flowcell,arch._meta.compress_suffix))
    dest_path_md5 = "{}.md5".format(dest_path)

    # If destination path already exists, check if it should be overwritten or, if md5sums match, it should be left as is and used
    if os.path.exists(dest_path):
        if not kw.get("force_overwrite",False):
            # Check md5sum
            if os.path.exists(dest_path_md5) and arch.app.cmd.verify_md5sum(dest_path_md5):
                arch.log.info("Run package already exists and md5sum matches, will not overwrite. Use --force-overwrite option to replace existing package")
                return dest_path
            elif os.path.exists(dest_path_md5):
                arch.log.warn("Run package already exists but md5sum does not match, will replace existing package")
            else:
                arch.log.warn("Run package already exists but no md5sum to compare against could be found, will replace existing package")
        else:
            arch.log.info("Run package already exists but --force-overwrite specified, will replace existing package")

    cmd = "tar {} --use-compress-program={} {}-cf {} -C {} {}".format(arch._meta.compress_opt,
                                                                arch._meta.compress_prog,
                                                                "--exclude-from={} ".format(excludes) if excludes else "",
                                                                dest_path,
                                                                root,
                                                                flowcell)
    # Run the compression
    arch.app.cmd.command(shlex.split(cmd), capture=True, ignore_error=False, cwd=workdir)

    # Calculate the md5sum
    arch.app.cmd.md5sum(dest_path)

    return dest_path

def upload_tarball(arch, tarball, remote_host=None, remote_path=None, remote_user=None, **kw):
    """Upload the tarball to the remote destination
    """
    if not remote_path:
        arch.log.error("A remote path must be specified in the config or on the command line")
        return False

    source_files = {'tarball': tarball,
                    'tarball_md5': "{}.md5".format(tarball)}

    arch.log.debug("Verifying that md5sum file {} exists".format(source_files['tarball_md5']))
    if not os.path.exists(source_files['tarball_md5']):
        arch.log.warn("md5 file {} does not exist".format(source_files['tarball_md5']))
        if not query_yes_no("Calculate md5 file and proceed?",
                            force=arch.pargs.force):
            return False

        # Calculate the md5sum
        arch.app.cmd.md5sum(source_files['tarball'])

    remote_location = "{}{}".format("{}@".format(remote_user) if remote_user else "",
                                    "{}:".format(remote_host) if remote_host else "")
    # Transfer the md5 file and tarball
    remote_files = {}
    for label in source_files.keys():
        remote_files[label] = "{}{}".format(remote_location,
                                            os.path.join(remote_path,os.path.basename(source_files[label])))
        arch.log.debug("Transferring {} to {}".format(source_files[label],remote_files[label]))
        arch.app.cmd.transfer_file(source_files[label],remote_files[label])

    # Verify the transfer on the remote side using fabric (if necessary)
    use_fabric = remote_host is not None and remote_host != "localhost"
    passed = False
    arch.log.debug("Verifying integrity of remote file {} after transfer".format(remote_files['tarball']))
    if use_fabric:
        # Verify the md5sum using fabric
        host, path = remote_files['tarball_md5'].split(':')
        result = execute(verify_upload,path,host=host)
        passed = result.get(host,False)
    else:
        passed = arch.app.cmd.verify_md5sum(remote_files['tarball_md5'])

    # If the verification was not successful, prompt to delete the corrupt files
    if not passed:
        arch.log.error("md5 sum of remote file {} does not match after transfer".format(remote_files['tarball']))
        if query_yes_no("Remove the corrupted remote file {}?".format(remote_files['tarball']),
                        force=arch.pargs.force):
            for path in remote_files.values():
                arch.log.info("removing {}".format(path))
                if use_fabric:
                    path = path.split(':')[-1]
                    execute(rm_file,path,host=host)
                else:
                    arch.app.cmd.safe_unlink(path)
        arch.log.error("Upload of {} to remote destination failed".format(source_files['tarball']))
    else:
        arch.log.info("{} uploaded to {} successfully".format(source_files['tarball'],remote_files['tarball']))

    if use_fabric:
        disconnect_all()

    return passed

def send_to_swestore(arch, tarball, swestore_path=None, remote_swestore=False, remote_user=None, remote_host=None, remote_path=None, **kw):
    """Method that sends a tarball to swestore using irods

    For now, the irods commands are just shell commands. These could probably be migrated to use the PyRods library.
    """

    passed = False
    if swestore_path is None:
        arch.log.error("A Swestore path must be specified in the config or with the --swestore-path command line option")
        return False

    # Execute on a remote host
    if remote_swestore:
        if not remote_path:
            arch.log.error("a remote path where the tarball is located must be specified in the config or on the command line")
            return False

        host = "{}{}".format("{}@".format(remote_user) if remote_user else "",
                             remote_host)
        remote_tarball = os.path.join(remote_path,os.path.basename(tarball))
        arch.log.info("executing the swestore archiving script on the remote host {}".format(host))
        result = execute(remote_swestore_archiving,remote_tarball,swestore_path,arch.pargs.clean,arch.pargs.dry_run,hosts=[host])
        passed = result.get(host,False)
        disconnect_all()

    # Run the script locally
    else:
        dir = os.path.dirname(tarball)
        fname = os.path.basename(tarball)

        command_kw = {'capture':True, 'ignore_error':False, 'shell':True}
        # Execute the archiving script
        try:
            cmd = "swestore_archive_run.sh {}{}{} {}".format("-d " if arch.pargs.clean_swestore else "",
                                                             "-n " if arch.pargs.dry_run else "",
                                                             tarball,
                                                             swestore_path)
            arch.log.info("sending {} to swestore {}".format(tarball,swestore_path))
            output = arch.app.cmd.command([cmd],**command_kw)
            arch.log.info(output)
            passed = True
        except:
            arch.log.error("archiving script failed")

    if not passed:
        arch.log.error("sending of {} to swestore  failed".format(tarball))
    else:
        arch.log.info("{} uploaded to {} in swestore successfully".format(tarball,swestore_path))

    return passed

@task
def verify_upload(remote_md5):
    """Verify the md5 sum of a remote file
    """
    # Go to the remote folder and execute the md5sum check
    remote_path = os.path.dirname(remote_md5)
    remote_fname = os.path.splitext(os.path.basename(remote_md5))[0]
    with settings(warn_only=True):
        with cd(remote_path):
            return run("md5sum -c {}".format(os.path.basename(remote_md5))).succeeded

@task
def rm_file(path):
    """Remove the supplied file
    """
    run("rm {}".format(path))

@task
def remote_swestore_archiving(tarball, swestore_path, clean, dry_run):
    """Run the swestore archiving code on the remote host
    """
    cmd = "pm archive swestore --tarball {} --send-to-swestore --swestore-path {}{}{}".format(tarball,
                                                                                              swestore_path,
                                                                                              " --clean" if clean else "",
                                                                                              " -n" if dry_run else "")
    with settings(warn_only=True):
        return run(cmd).succeeded


