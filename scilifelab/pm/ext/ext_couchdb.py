"""Couchdb extension"""

import os
import sys
import couchdb

from cement.core import backend, handler, hook

from scilifelab.pm.core import command
from scilifelab.utils.http import check_url
from scilifelab.utils.timestamp import utc_time

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

        url = None
        """The database url"""

        user = None
        """The user"""
        
        pw = None
        """The password used to connect"""

        views = {}
        """Temporary views for speeding up update functions"""

    ### Must be present
    def command(self):
        pass

    def connect(self, url, port="5984"):
        def runpipe():
            self._meta.url="http://{}:{}".format(url,port)
            if not check_url(self._meta.url):
                self.app.log.warn("Connecting to server at {} failed. No such url." % self._meta.url)
                return
            self._meta.conn = couchdb.Server(url=self._meta.url)
            self.app.log.info("Connecting to server at {} succeeded".format(self._meta.url))
        return self.dry("Connecting to database @{}:{}".format(url, port), runpipe)

    def db(self, dbname):
        """Get database <dbname>.

        :param dbname: Database name

        :returns: database on success, otherwise False
        """
        def runpipe():
            try:
                db = self._meta.conn[dbname]
            except:
                self.app.log.warn("No such database {} ".format(dbname))
                return False
            return db
        return self.dry("Retrieving database {} from {}".format(dbname, self._meta.url), runpipe)

    def save(self, dbname, obj, update_fn=None):
        """Save/update database object <obj> in database <dbname>. If
        <obj> already exists and <update_fn> is passed, update will
        only take place if object has been modified

        :param dbname: database name
        :param obj: database object to save
        :param update_fn: function that operates on object and makes sure it doesn't already exist
        """
        def runpipe():
            db = self.db(dbname)
            if not update_fn:
                db.save(obj)
            else:
                new_obj = update_fn(db, obj)
                if not new_obj is None:
                    self.app.log.info("Saving object {} with id {}".format(repr(new_obj), new_obj["_id"]))
                    db.save(new_obj)
                else:
                    self.app.log.info("Object {} present and not in need of updating".format(repr(obj)))
        return self.dry("Saving object {}".format(repr(obj)), runpipe)

    def get_view(self, dbname, design, name):
        """Get view from a database <dbname> with design document <design>, named <name>

        :param dbname: database name
        :param design: design document
        :param name: view name
        """
        db = self.db(dbname)
        return db.view("{}/{}".format(design, name))


    def _view(self, dbname, key):
        """Create or retrieve a 'view' in database <dbname> that maps <field> to _id and _rev.
        
        :param dbname: database name
        :param key: database key

        """
        db = self.db(dbname)
        k = "{}_{}".format(dbname, key)
        if not self._meta.views.has_key(k):
            self.app.log.info("generating view in database '{}' for key '{}'".format(dbname, key))
            self._meta.views[k] = {}
            for dbid in db:
                dbobj = db.get(dbid)
                value = dbobj.get(k)
                self._meta.views[k][value] = (dbobj.get("_id"), dbobj.get("_rev"))
        return self._meta.views[k]

    def _get(self, dbname, key, query):
        """Get a field from a view where field == <query>.

        :param dbname: database name
        :param key: database key
        :param query: query
        """
        db = self.db(dbname)
        view = self._view(dbname, key)
        return view.get(query, None)
                
def add_shared_couchdb_options(app):
    """
    Adds shared couchdb arguments to the argument object.
    
    :param app: The application object.
    
    """
    app.args.add_argument('--url', help="Database url (excluding http://)", nargs="?", type=str)
    app.args.add_argument('--port', help="Database port. Default 5984", nargs="?", default="5984", type=str)
    app.args.add_argument('--dbname', help="Database name", default=None, type=str)
    app.args.add_argument('--user', help="Database user", nargs="?", default=None, type=str)
    app.args.add_argument('--password', help="Database password", default=None, type=str)


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
