#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

setup(name = "scilifelab",
      version = "0.1",
      author = "SciLife",
      author_email = "genomics@scilifelab.se",
      description = "Useful scripts for use at SciLife",
      license = "MIT",
      scripts = ['scripts/runsizes.py',
                 'scripts/bcbb_helpers/run_bcbb_pipeline.py',
                 'scripts/bcbb_helpers/report_to_gdocs.py',
                 'scripts/bcbb_helpers/process_run_info.py',
                 'scripts/pm'],
      install_requires = [
          "bcbio-nextgen >= 0.2",
          "drmaa >= 0.5",
	  "sphinx >= 1.1.3",
	  "couchdb >= 0.8",
          "reportlab >= 2.5",
          "cement >= 2.0.2",
      ],
      test_suite = 'nose.collector',
      #packages=find_packages(exclude=['tests']),
      packages=['scilife'],
      #package_dir = {'pmtools':'scilife/project_management/pmtools'},
      #package_data = {'pmtools':['scilife/project_management/pmtools/templates/tpl/make/*']}
      )
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
