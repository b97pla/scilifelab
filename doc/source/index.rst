.. scilifelab documentation master file, created by
   sphinx-quickstart on Fri Aug 31 10:26:17 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SciLifeLab modules
===================================

SciLifeLab modules is a collection of python modules, mainly used for
maintenance, production, and analysis of next-generation sequencing
data. 


Project management tools (pm)
-----------------------------

*pm* is a client application for performing common project
management tasks of next-generation sequencing data. It is primarily
intented for use at `Science for Life Laboratory.
<http://www.scilifelab.se>`_ The interface is based on `Cement.
<http://builtoncement.com/2.1/>`_ 

Features include:

 * Submission of jobs via DRMAA
 * Compression of project files
 * Cleaning of analysis directories
 * Running pre-defined workflows that work as add-ons to the bcbio
   pipeline


Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation

API
---

.. toctree::
   :maxdepth: 1

   api/index
   pm/index
