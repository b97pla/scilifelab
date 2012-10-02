Configuration
=============

In order to run *pm*, there must be a configuration file located at
"~/.pm/pm.conf". Since this module is geared towards use at SciLife,
one might wonder why there is a config file in the first place since
we know where everything is located. Two reasons come to mind:

1. locations may actually change
2. locations shouldn't be hardcoded anyway

Point 2 does require some care though, since pointing the config paths
to wrong locations may cause undesirable effects.

An example configuration file
-----------------------------

.. code-block:: text

   [archive]
   root = /home/username/data/archive

   [production]
   root = /home/username/data/production

   [project]
   root = /home/username/data/projects
   repos = /home/username/data/repos
   finished = /home/username/data/projects/finished

   [log]
   file=/home/user/pm.log
