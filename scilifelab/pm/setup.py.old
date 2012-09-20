#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

setup(name = "pm",
      version = "0.1.0",
      author = "Per Unneberg",
      author_email = "per.unneberg@scilifelab.se",
      description = "Project management tools",
      namespace_packages=["pmtools"],
      license = "MIT",
      scripts = ['scripts/pm'],
      install_requires = [
          "drmaa >= 0.5",
	  "sphinx >= 1.1.3",
	  "couchdb >= 0.8",
          "cement >= 2.0.2",
      ],
      test_suite = 'nose.collector',
      packages=find_packages(),
      package_dir = {'pmtools':'pmtools'},
      package_data = {'pmtools':['templates/tpl/make/*']}
      )
#os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
