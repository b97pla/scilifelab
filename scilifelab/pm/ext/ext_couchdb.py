"""Couchdb extension."""
from cement.core import hook

def add_shared_couchdb_options(app):
    """
    Adds shared couchdb arguments to the argument object.
    
    :param app: The application object.
    
    """
    user = None
    url = None
    password = None
    if app.config.has_option("db", "user"):
        user = app.config.get("db", "user") 
    if app.config.has_option("db", "password"):
        password = app.config.get("db", "password") 
    if app.config.has_option("db", "url"):
        url = app.config.get("db", "url") 
    group = app.args.add_argument_group('couchdb', 'Options for couchdb connections')
    group.add_argument('--url', help="Database url (excluding http://). Default '{}'".format(url), default=url, nargs="?", type=str)
    group.add_argument('--port', help="Database port. Default 5984", nargs="?", default="5984", type=str)
    group.add_argument('--username', help="Database user. Default '{}'".format(user), nargs="?", default=user, type=str)
    group.add_argument('--password', help="Database password.", default=password, type=str)

def load():
    """Called by the framework when the extension is 'loaded'."""
    hook.register('post_setup', add_shared_couchdb_options)
