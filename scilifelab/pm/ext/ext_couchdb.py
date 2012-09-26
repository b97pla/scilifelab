"""Couchdb extension"""

import os
import sys
import couchdb

from cement.core import backend, handler, hook

from scilifelab.pm.core import command
from scilifelab.utils.http import check_url

LOG = backend.minimal_logger(__name__)

class CouchdbCommandHandler(command.CommandHandler):
    """ 
    This class is an implementation of the :ref:`ICommand
    <scilifelab.pm.core.command>` interface.
    """    

    class Meta:
        """Handler meta-data"""
        
        interface = command.ICommand
        """The interface that this class implements."""

        label = 'couchdb'
        """The string identifier of this handler."""

        conn = None
        """The database connection"""

    def connect(self, url, port="5984"):
        def runpipe():
            db_url="http://{}:{}".format(url,port)
            if not check_url(db_url):
                self.app.log.warn("Connecting to server at {} failed. No such url." % db_url)
                return
            self.conn = couchdb.Server(url=db_url)
            self.app.log.info("Connecting to server at {} succeeded".format(url))
        return self.dry("Connecting to database @{}:{}".format(url, port), runpipe)

    def db(self, dbname=None):
        if not dbname:
            return
        try:
            db = self.conn[dbname]
        except:
            self.app.log.warn("No such database {}".format(dbname))
            return False
        return db
        

    def command(self, cmd_args, capture=True, ignore_error=False, cwd=None, **kw):
        def runpipe():
            pass
        return self.dry("Connecting to database @{}:{}".format(url, port), runpipe)

def add_shared_couchdb_options(app):
    """
    Adds shared couchdb arguments to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('url', help="Database url (excluding http://)", nargs="?", type=str)
    app.args.add_argument('--port', help="Database port. Default 5984", nargs="?", default="5984", type=str)
    app.args.add_argument('--dbname', help="Database name", default=None, type=str)

def add_couchdb_option(app):
    """
    Adds the '--couchdb' argument to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--couchdb', dest='output_handler', 
                          action='store_const', help='toggle couchdb output handler', const='couchdb')


def set_couchdb_handler(app):
    """
    Overrides the configured command handler if ``--couchdb`` is passed at the
    command line.
    
    :param app: The application object.
    
    """
    if '--couchdb' in app._meta.argv:
        app._meta.cmd_handler = 'couchdb'
        app._setup_cmd_handler()

def load():
    """Called by the framework when the extension is 'loaded'."""
    hook.register('post_setup', add_couchdb_option)
    hook.register('post_setup', add_shared_couchdb_options)
    hook.register('pre_run', set_couchdb_handler)
    handler.register(CouchdbCommandHandler)
