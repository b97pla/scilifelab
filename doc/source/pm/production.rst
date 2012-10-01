Production module commands
==========================

Provide functionality for project management.

Commands:
^^^^^^^^^

.. code:: text      

       ls
         list contents
       init
         initialize a project folder
       add
         add boilerplate code
       compress
         compress files
       clean      
         remove files
       du          
         calculate disk usage
       deliver     
         deliver project data to customer

Synopsis
---------

The following command creates a directory in the project root
named j_doe_00_00. The '-g' flag adds a git directory to the
project repos, and initializes the project subdirectory j_doe_00_00_git 
for use with git.

.. code:: bash

   pm project init j_doe_00_00 -g

FIXME: Boilerplate code can be added to the project by running

.. code:: bash

   pm project add j_doe_00_00

The boilerplate code includes makefiles, sbatch templates, and documentation
templates.
