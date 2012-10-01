#!/usr/bin/env python
from setuptools import setup, find_packages
import sys
import os
import glob

setup(name = "scilifelab",
      version = "0.2",
      author = "Science for Life Laboratory",
      author_email = "genomics_support@scilifelab.se",
      description = "Useful scripts for use at SciLifeLab",
      license = "MIT",
      namespace_packages=["scilifelab"],
      scripts = glob.glob('scripts/*.py') + ['scripts/pm', 'scripts/pm-deliver'],
                 #'scripts/bcbb_helpers/*.py,
      install_requires = [
        "bcbio-nextgen >= 0.2",
        "drmaa >= 0.5",
        "sphinx >= 1.1.3",
        "couchdb >= 0.8",
        "reportlab >= 2.5",
        "cement >= 2.0.2",
        "mock",
        "PIL"
        ],
      test_suite = 'nose.collector',
      packages=['scilifelab'],
      ## package_data: install data/templates needed by modules
      package_data = {'scilifelab':[
            'pm/templates/tpl/make/*',
            'data/grf/*',
            #'sbatch/*' # should be deprecated by python-drmaa, not even installing it...
            ]}
      )

os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
