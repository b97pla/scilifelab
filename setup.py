#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

<<<<<<< HEAD
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
	  "couchdb >= 0.8",
          "reportlab >= 2.5",
      ])
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
=======
version = '0.0'

setup(name='scilifelab',
      version=version,
      description="SciLifeLab life sciences facility platform python scripts.",
      long_description="""Science For Life Laboratory python scripts""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='scilifelab science',
      author='Science for Life Lab staff',
      author_email='genomics@scilifelab.se',
      url='http://scilifelab.se',
      license='GPLv3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          "cement"
      ],
      entry_points={
        'console_scripts': [
            'pm = scilifelab.bin.pm'
         ]
      }
   )
>>>>>>> 7332c286b0a31129ee9aa2a34760e4b81352228e
