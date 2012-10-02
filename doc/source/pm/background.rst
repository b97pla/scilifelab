Background
----------

The *pm* client application grew out of an increasing frustration over
the lack of standardized tools for our production environment and (in
particular) our best practice environment. Sure, there are a number of
scripts that do this and that, but it's not always easy to know where
they are, or what they do. There was a need for a common entry point
that gathers all relevant functionality for standardized, repeatable
tasks, such as compressing/cleaning files in production folders,
setting up project folders with a common directory structure, or
listing the contents of a project in a given way. In addition, any
operation that tampers with the files in the production environment
should be logged.

*pm* is an attempt at such a common interface. It is built upon the
python module *cement*, an advanced CLI Application Framework. Apart
from providing an easy way to adding command line applications, it
also has a nice logging facility. It is also possible to create
plugings and extensions, making it very easy to add customized code.


