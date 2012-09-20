import os
from cement.core import backend, handler, output
from cement.utils import test
from scilifelab.pm import PmApp

## Set default configuration
filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
config_defaults = backend.defaults('analysis', 'archive', 'config', 'project','log')
config_defaults['analysis']['root']  = os.path.join(filedir, "data", "analysis")
config_defaults['archive']['root']  = os.path.join(filedir, "data", "archive")
config_defaults['project']['root']  = os.path.join(filedir, "data", "projects")
config_defaults['project']['repos']  = os.path.join(filedir, "data", "repos")
config_defaults['config']['ignore'] = ["slurm*", "tmp*"]
config_defaults['log']['level']  = "INFO"
config_defaults['log']['file']  = os.path.join(filedir, "data", "log", "pm.log")


def safe_makedir(dname):
    """Make directory"""
    if not os.path.exists(dname):
        try:
            os.makedirs(dname)
        except OSError:
            if not os.path.isdir(dname):
                raise
    else:
        print "Directory %s already exists" % dname
    return dname


## Output handler for tests
class PmTestOutputHandler(output.CementOutputHandler):
    class Meta:
        label = 'pmtest'

    def render(self, data, template = None):
        for key in data:
            if data[key]:
                print "{} => {}".format(key, data[key].getvalue())

## Testing app
class PmTestApp(PmApp):
    class Meta:
        argv = []
        config_files = []
        config_defaults = config_defaults
        output_handler = PmTestOutputHandler

    # def setup(self):
    #     super(PmTestApp, self).setup()
    #     self._output_data = dict(stdout=[], stderr=[])


## Main pm test 
class PmTest(test.CementTestCase):
    app_class = PmTestApp
    app = None
    OUTPUT_FILES = []

    def setUp(self):
        self._clean()

    def _clean(self):
        pass

    def _run_app(self):
        try:
            self.app.setup()
            self.app.run()
            self.app.render(self.app._output_data)
        finally:
            self.app.close()
            
        
