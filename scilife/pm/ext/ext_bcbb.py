"""
bcbb extension
==============


"""
import os
import sys

from cement.core import backend, handler, hook

from scilife.pm.core import naming

LOG = backend.minimal_logger(__name__)

class BcbbNamingHandler(naming.NamingHandler):
    """
    This class is an implementation of the :ref:`ICommand
    <pmtools.core.command>` interface.
    """    

    class Meta:
        """Handler meta-data"""
        
        interface = naming.INaming
        """The interface that this class implements."""

        label = 'bcbb'
        """The string identifier of this handler."""


    ext_format = dict(
        nophix_fastq_txt = '{}_{}_{}_{}_nophix_{}_'
        )
