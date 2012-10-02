"""Pm Output Handler"""
import sys

from cement.core import output

class PmOutputHandler(output.CementOutputHandler):
    """
    Main Pm output handler.

    """
    class Meta:
        label = 'pmout'

    def render(self, data, template = None):
        """
        Render output data stored in cStringIO objects in data.

        :param data: dictionary with keys <stdout> and <stderr>
        :param template: template output. Currently not implemented.
        """
        if data["stdout"].getvalue():
            print >> sys.stdout, data["stdout"].getvalue()
        if data["stderr"].getvalue():
            print >> sys.stderr, data["stderr"].getvalue()
