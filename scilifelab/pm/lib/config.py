"""Configuration settings"""

import os
from cement.core import backend

config_defaults = backend.defaults('production', 'archive', 'config', 'project','log', 'distributed')
config_defaults['production']['root']  = None
config_defaults['archive']['root']  = None
config_defaults['project']['root']  = None
config_defaults['project']['repos']  = None
config_defaults['project']['finished']  = None
config_defaults['config']['ignore'] = ["slurm*", "tmp*"]
config_defaults['log']['level']  = "INFO"
config_defaults['log']['file']  = os.path.join(os.getenv("HOME"), "log", "pm.log")
config_defaults['distributed']['jobaccount'] = None
