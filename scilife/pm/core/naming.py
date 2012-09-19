"""Naming interface"""

import os
import sys

from cement.core import interface, handler

def naming_interface_validator(cls, obj):
    members = [
        '_setup',
        'ext_format',
        ]

class INaming(interface.Interface):

    class IMeta:
        """Interface meta-data"""

        label = 'naming'
        """Naming interface"""

        validator = naming_interface_validator
        """Interface validator function"""
        
    Meta = interface.Attribute('Handler meta-data')
    ext_format = interface.Attribute('Extension to format mapping')

    def _setup(app_obj):
        """
        The _setup function is called during application initialization and
        must 'setup' the handler object making it ready for the framework
        or the application to make further calls to it.
        
        :param app_obj: The application object. 
                                
        """

class NamingHandler(handler.CementBaseHandler):
    """
    Base class that all Naming Handlers should sub-class from.
    """
    class Meta:
        """
        Handler meta-data
        """

        label = None
        """String identifier."""

        interface = INaming
        """The interface"""

    def __init__(self, *args, **kw):
        super(NamingHandler, self).__init__(*args, **kw)


