"""classes for full test"""
import os
from cement.core import backend, handler, output
from cement.utils import test
from scilifelab.pm import PmApp

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

config_defaults = backend.defaults('production', 'archive', 'config', 'project','log', 'runqc', 'db')
config_defaults['production']['root']  = os.path.join(filedir, "data", "production")
config_defaults['archive']['root']  = os.path.join(filedir, "data", "archive")
config_defaults['project']['root']  = os.path.join(filedir, "data", "projects")
config_defaults['project']['repos']  = os.path.join(filedir, "data", "repos")
config_defaults['runqc']['root']  = os.path.join(filedir, "data", "archive")
config_defaults['runqc']['production']  = os.path.join(filedir, "data", "production")
config_defaults['runqc']['rsync_sample_opts'] = "-amnv  --include=*/ --include='[0-9][0-9]_\*.txt' --include='bcbb_software_versions.txt' --include='*.yaml*' --include='*screen.txt' --include='*summary.txt' --include='*fastqc_data.txt' --include='*metrics' --exclude='tmp' --exclude='*'"
config_defaults['config']['ignore'] = ["slurm*", "tmp*"]
config_defaults['log']['level']  = "INFO"
config_defaults['log']['file']  = os.path.join(filedir, "data", "log", "pm.log")
config_defaults['db']['url'] = "localhost"
config_defaults['db']['user'] = "u"
config_defaults['db']['password'] = "p"
config_defaults['db']['samples'] = "samples-test"
config_defaults['db']['flowcells'] = "flowcells-test"
config_defaults['db']['projects'] = "projects-test"

## Output handler for tests
class PmTestOutputHandler(output.CementOutputHandler):
    class Meta:
        label = 'pmfulltest'

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

class PmFullTest(test.CementTestCase):
    app_class = PmTestApp
    app = None

    def setUp(self):
        pass

    def _run_app(self):
        try:
            self.app.setup()
            with self.app.log.log_setup.applicationbound():
                self.app.run()
                self.app.render(self.app._output_data)
        finally:
            self.app.close()
        
