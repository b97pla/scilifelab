#!/usr/bin/env python
"""Setup file and install script SciLife python scripts.
"""
from setuptools import setup, find_packages
import os

setup(name = "scilifelab",
      version = "0.1",
      author = "SciLife",
      author_email = "genomics@scilifelab.se",
      description = "Useful scripts for use at SciLife",
      license = "MIT",
      scripts = ['scripts/project_management.py',
                 'scripts/runsizes.py',
                 'scripts/bcbb_helpers/run_bcbb_pipeline.py',
                 'scripts/bcbb_helpers/process_run_info.py'],
      install_requires = [
          "bcbio-nextgen >= 0.2",
          "drmaa >= 0.5",
	  "sphinx >= 1.1.3",
      ])
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
