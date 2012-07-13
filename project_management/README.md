Title: pm documentation
Author: Per Unneberg <per.unneberg@scilifelab.se>
Date: 2012-07-13 08:17:12 CEST

pm documentation 
`========`

# Synopsis #

Project management tools - wrappers and help functions to facilitate
working with projects at scilife

`
pm
pm archive ls
pm archive runinfo FLOWCELLID -t
pm project ls -h
pm project ls PROJECTID
pm analysis status
pm analysis status FLOWCELLID
pm analysis clean FLOWCELLID -p projectid
`

# Motivation #

I have had a growing frustration over the lack of standardized tools
for our production environment and (in particular) our best practice
environment. Sure, there are a number of scripts that do this and
that, but it's not always easy to know where they are, or what they
do. I felt we needed a common entry point that gathers all relevant
functionality for standardized, repeatable tasks, such as
compressing/cleaning files in analysis folders, setting up project
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

# Modules #

## archive ##

## analysis ##

## project ##


# References #

[^1]: http://cement.readthedocs.org
