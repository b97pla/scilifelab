Installation
------------

The package is shipped with a *pavement.py* file that is used to
generate a *setup.py* file used for installation. To install, simply
change directory to the root directory ('project_management') and
issue the following commands:

.. code-block:: shell

   paver generate_setup
   python setup.py install

To make documentation files, cd to the 'doc' subdirectory and type 

.. code-block:: shell

   make html
