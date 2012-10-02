Scilifelab modules
==================

Installation
------------

Installation is as simple as

.. code:: bash

   python setup.py install

If you are running several virtual environments, where one (e.g.
`devel`) is used for development, you can install a development
version by running

.. code:: bash

   workon devel
   python setup.py develop


Documentation
--------------

Docs are located in the `doc` directory. To install, cd to `doc` and
run

.. code:: bash

   make html

Documentation output is found in `build`.


Running the tests
-----------------

The modules are shipped with a number of unit tests located in the
`tests` directory. To run a test, issue the command

.. code:: bash

   python setup.py nosetests

or if you want to run individual tests, cd to `tests` and run (for
example)

.. code:: bash

   nosetests -v -s test_db.py

