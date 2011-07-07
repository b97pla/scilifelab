#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts at SciLife Lab.
"""
from setuptools import setup, find_packages

setup(name = "scilife-ngs",
      version = "0.1",
      author = "Per Unneberg",
      author_email = "per.unneberg@scilifelab.se",
      description = "Next generation sequencing analysis workflows",
      license = "BSD",
      url = "",
      namespace_packages = ["scilife"],
      packages = find_packages(),
      scripts = ['scripts/exome_pipeline.py',
                 ],
      package_data = {
          'config' : ['*.yaml'],
          },
      install_requires = [
          "biopython >= 1.56",
          "Mako >= 0.3.6",
          "PyYAML >= 3.09",
          "Logbook >= 0.3",
          "pysam >= 0.4.1",
          "bcbio >= 0.1"
          #"fabric >= 1.0.1",
          #"celery >= 2.2.7",
          #"rpy2 >= 2.0.7"
      ])
