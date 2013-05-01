Adding code
============

The following sections show brief examples of how to add commands
making use of classes defined in *pm*. See [^1] for more detailed
information.

Controllers
------------

.. note:: adding a controller requires modification to the main
          application script *pm*.

Assume you want to create a controller *mycommand*. Create a file
called *mycommand.py* and place it in *pmtools/controllers*. Create a
doc string, and import the follow two modules:

.. code-block:: python

    from cement.core import controller
    from pmtools import AbstractBaseController

The *AbstractBaseController* is an interface, ensuring that all *pm*
controllers behave similarly. The minimum boilerplate code needed to
define the controller is


.. code-block:: python

    ## Main mycommand controller
    class MycommandController(AbstractBaseController):
        """
        Functionality for mycommand.
        """
        class Meta:
            label = 'mycommand'
            description = 'Manage mycommand'
	    stacked_on = None

        @controller.expose(hide=True)
        def default(self):
            pass

To add subcommands, add functions decorated with @controller.expose:

.. code-block:: python

    @controller.expose(help="My subcommand")
    def mysubcommand(self):
	print "Mysubcommand"
		
That's all there is to it.

Before running, we need to modify the main application script. First
you need to import the newly defined command:

.. code-block:: python

    from pmtools.controller.mycommand import MycommandController
	
Then, before the application is setup (`app.setup()`) the command
needs to be registered:

.. code-block:: python

    handler.register(MycommandController)

The new command can now be accessed as 

.. code-block:: bash

    pm mycommand
    pm mycommand mysubcommand


Subcontrollers
--------------

The main controllers are unstacked, i.e. their arguments are specific
to each controller. However, one can also add stacked controllers that
add arguments to the main controllers.

There also is an interface to subcontrollers called *SubController*.
To add a subcontroller to *mycommand.py*, add

.. code-block:: python

    class Mysubcommand2Controller(SubController):
        class Meta:
	    label = 'mysubcommand2-ctrl'
	    stacked_on = 'mycommand'
	    description = 'Mysubcommand2 controller'
	    interface = controller.IController
	    arguments = [
	    	(['-f', '--foo'], dict(help="foo argument", default=False, action="store_true"))
			]

    @controller.expose(help="Mysubcommand2 help")
    def mysubcommand2(self):
        print "mysubcommand2"
		


Plugins
-------


