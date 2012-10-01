Title: pm documentation

Author: Per Unneberg <per.unneberg@scilifelab.se>

Date: 2012-07-13 08:17:12 CEST

# pm documentation #

Project management tools - wrappers and help functions to facilitate
working with projects at scilife

## Synopsis ##

	pm
    pm archive ls
    pm archive runinfo FLOWCELLID -t
    pm project ls -h
    pm project ls PROJECTID
    pm production status
    pm production status FLOWCELLID
    pm production clean FLOWCELLID -p projectid
    
# Motivation #

I have had a growing frustration over the lack of standardized tools
for our production environment and (in particular) our best practice
environment. Sure, there are a number of scripts that do this and
that, but it's not always easy to know where they are, or what they
do. I felt we needed a common entry point that gathers all relevant
functionality for standardized, repeatable tasks, such as
compressing/cleaning files in production folders, setting up project
folders with a common directory structure, or listing the contents of
a project in a given way. In addition, any operation that tampers with
the files in the production environment should be logged.

**pm** is my attempt at a common interface. It is built upon the
python module **cement**[^1], an advanced CLI Application Framework.
Apart from providing an easy way to adding command line applications,
it also has a nice logging facility, and customizable. It is also
possible to create plugings, a feature which I have not yet tried out.

## Why the config file? ##

Since this module is geared towards use at SciLife, one might wonder
why there is a config file in the first place since we know where
everything is located. Two reasons come to mind:

1. locations may actually change
2. locations shouldn't be hardcoded anyway

Point 2 does require some care though, since pointing the config paths
to wrong locations may cause undesirable effects.

# Installation #

First, clone the repository and cd to the install location. In the
root directory you should see the file *pavement.py*. Run

`paver generate_setup`

which will generate *setup.py*. Then type

`python setup.py install`

and you're done.

# Adding code #

The following sections show brief examples of how to add commands
making use of classes defined in *pmtools*. See [^1] for more detailed
information.

## Controllers ##

NOTE: adding a controller requires modification to the main
application script *pm*.

Assume you want to create a controller *mycommand*. Create a file
called *mycommand.py* and place it in *pmtools/controllers*. Create a
doc string, and import the follow two modules:

    from cement.core import controller
    from pmtools import AbstractBaseController

The *AbstractBaseController* is an interface, ensuring that all *pm*
controllers behave similarly. The minimum boilerplate code needed to
define the controller is

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

    @controller.expose(help="My subcommand")
    def mysubcommand(self):
		print "Mysubcommand"
		
That's all there is to it.

Before running, we need to modify the main application script. First
you need to import the newly defined command:

    from pmtools.controller.mycommand import MycommandController
	
Then, before the application is setup (`app.setup()`) the command
needs to be registered:

    handler.register(MycommandController)

The new command can now be accessed as 

    pm mycommand
	pm mycommand mysubcommand

### Subcontrollers ###

The main controllers are unstacked, i.e. their arguments are specific
to each controller. However, one can also add stacked controllers that
add arguments to the main controllers.

There also is an interface to subcontrollers called *SubController*.
To add a subcontroller to *mycommand.py*, add

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
		


## Plugins ##

TODO: I think maybe the above section should be restricted to a
well-defined set of commands, and that we rather work with plugins to
extend functionality

# References #

[^1]: http://cement.readthedocs.org
