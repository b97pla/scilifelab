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
      scripts = ['scripts/project_management.py'],
      install_requires = [
          "bcbio-nextgen >= 0.2",
      ])
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
