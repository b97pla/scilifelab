"""scilifelab.utils.halo module"""
import os
from cStringIO import StringIO
from mako.template import Template
from scilifelab.experiment.project import find_samples
import scilifelab.log
from scilifelab.utils.dry import dry_write

LOG = scilifelab.log.minimal_logger(__name__)
FILEPATH=os.path.dirname(os.path.realpath(__file__))

# FIXME: change when moving from experimental
HALORC=os.path.join(FILEPATH, os.pardir, os.pardir, "experimental", "halorc")
HALOSCRIPT=os.path.join(FILEPATH, os.pardir, os.pardir, "experimental", "halo_pipeline.sh")

TEMPLATEPATH = os.path.join(FILEPATH, os.pardir, "data", "templates", "halo")
PROJECTTEMPLATE = Template(filename=os.path.join(TEMPLATEPATH, "projectrc.mako"))

def run_halo(path=None, project=None, batch_size=8, **kw):
    """Run halo application. Setup parameter files and call
    halo_pipeline.sh script.

    :param project: project name
    :param batch_size: number of samples to run in each project config file
    """
    plist = sorted(find_samples(path, **kw))
    plist_chunks=[plist[x:x+batch_size] for x in xrange(0, len(plist), batch_size)]
    i = 0
    param_list = []
    for pl in plist_chunks:
        i += 1
        param = {'cl':None, 'platform_args':None, 'workingDirectory':None}
        d = {'samples' : '"{}"'.format(" ".join([os.path.basename(x) for x in pl])),
             'indir' : path,
             'baits_file' : kw.get('baits', ""),
             'targets_file' : kw.get('targets', ""),
             'target_region' : kw.get('target_region', "")
             }
        outfile = os.path.join(path, "{}_{}_halo.projectrc".format(project, i))
        if kw.get("setup", False):
            dry_write(outfile, PROJECTTEMPLATE.render(**d), dry_run=kw.get("dry_run", False))
        if not os.path.exists(outfile):
            LOG.warn("No such configuration file {}; rerun command with '--setup' option")
            return []
        if kw.get("config", None) and os.path.basename(outfile) != kw.get("config", None):
            continue
        param['cl'] = [HALOSCRIPT, "-c", HALORC, outfile]
        label = '{}_halo_{}'.format(project[0:3].replace(".", "_"), i)
        param['platform_args'] = ['--output', os.path.join("{}.out".format(label)),
                                  '--error', os.path.join("{}.err".format(label)),
                                  '--job-name', label]
        param['workingDirectory'] = os.path.dirname(outfile)
        param_list.append(param)
    return param_list
