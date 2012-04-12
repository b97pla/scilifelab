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
      scripts = ['scripts/analyze_quality_recal.py',
                           ],
      package_data = {
          'config' : ['*.yaml'],
          },
      install_requires = [
          "bcbb >= 0.3a",
      ])
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
