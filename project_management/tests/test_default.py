import os
from cement.core import backend, handler
from cement.utils import test
from pmtools import PmApp

## Set default configuration
config_defaults = backend.defaults('analysis', 'archive', 'config', 'project','log')
config_defaults['analysis']['root']  = os.path.join(os.path.abspath(os.getcwd()), "data", "analysis")
config_defaults['archive']['root']  = os.path.join(os.path.abspath(os.getcwd()), "data", "archive")
config_defaults['project']['root']  = os.path.join(os.path.abspath(os.getcwd()), "data", "projects")
config_defaults['project']['repos']  = os.path.join(os.path.abspath(os.getcwd()), "data", "repos")
config_defaults['config']['ignore'] = ["slurm*", "tmp*"]
config_defaults['log']['level']  = "INFO"
config_defaults['log']['file']  = os.path.join(os.path.abspath(os.getcwd()), "data", "log", "pm.log")

## Testing app
class PmTestApp(PmApp):
    class Meta:
        argv = []
        config_files = []
        config_defaults = config_defaults


## Main pm test 
class PmTest(test.CementTestCase):
    app_class = PmTestApp
    app = None
    OUTPUT_FILES = []

    def setUp(self):
        self._clean()

    def _clean(self):
        print "cleaning"

    def _run_app(self):
        try:
            self.app.setup()
            self.app.run()
        finally:
            self.app.close()
        
