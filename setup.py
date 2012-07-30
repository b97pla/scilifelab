#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

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
