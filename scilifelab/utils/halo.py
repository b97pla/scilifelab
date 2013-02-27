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
    output_data = {'stdout':StringIO(), 'stderr':StringIO()}    
    plist = sorted(find_samples(path, **kw))
    plist_chunks=[plist[x:x+batch_size] for x in xrange(0, len(plist), batch_size)]
    i = 1
    for pl in plist_chunks:
        d = {'samples' : '"{}"'.format(" ".join([os.path.basename(x) for x in pl])),
             'indir' : path,
             'baits_file' : kw.get('baits', ""),
             'targets_file' : kw.get('targets', ""),
             'target_region' : kw.get('target_region', "")
             }
        outfile = os.path.join(path, "{}_{}_halo.projectrc".format(project, i))
        i += 1
        dry_write(outfile, PROJECTTEMPLATE.render(**d), dry_run=kw.get("dry_run", False))
        cl = [HALOSCRIPT, "-c", HALORC, outfile]
        print cl
    
    return output_data    
