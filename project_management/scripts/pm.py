import os
from cement.core import backend, foundation, hook, handler
import argparse
from pmtools import PmController
from pmtools.controller.project import ProjectController
from pmtools.controller.clean import CleanController
from pmtools.controller.compress import CompressController
from pmtools.controller.ls import LsController, RunInfoController

CONFIGFILE=os.path.join(os.getenv("HOME"), ".pm.conf")

## Tried to set this in controller subclass but doesn't work?!?
## FIXME: move to separate config.py module
defaults = backend.defaults('analysis', 'archive', 'config' 'project','log')
defaults['analysis']['root']  = None
defaults['archive']['root']  = None
defaults['project']['root']  = None
defaults['projects']['repos']  = None

defaults['config']['ignore'] = ["slurm*", "tmp*"]
defaults['log']['level']  = "INFO"
defaults['log']['file']  = os.path.join(os.getenv("HOME"), "log", "pm.log")


app = foundation.CementApp('pm', base_controller=PmController, config_defaults = defaults, extensions = ['json'])
    
try:
    handler.register(LsController)
    handler.register(RunInfoController)
    handler.register(CleanController)
    handler.register(CompressController)
    app.setup()

    try:
        app.config.parse_file(CONFIGFILE)
    except:
        app.log.warn("No config file %s; please create and set relevant config sections" % CONFIGFILE)
        sys.exit()
    app.run()

finally:
    app.close()
