"""Database module"""
import os
import sys
import couchdb

from scilifelab.log import minimal_logger
from scilifelab.utils.http import check_url

class ConnectionError(Exception):
    """Exception raised for connection errors.


    :param msg: the error message
    """

    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg

    def __str__(self):
        return self.msg

class Database(object):
    """Main database connection object for noSQL databases"""

    def __init__(self, **kwargs):
        self.con = None
        self.log = minimal_logger(repr(self))
        self.connect(**kwargs)

    def connect(self, **kwargs):
        pass

    def save(self, **kwargs):
        pass

    def __repr__(self):
        return "{}".format(self.__class__)

## From http://stackoverflow.com/questions/8780168/how-to-begin-writing-a-python-wrapper-around-another-wrapper
class Couch(Database):
    _doc_type = None
    _update_fn = None

    def __init__(self, log=None, url="localhost", **kwargs):
        self.db = None
        self.url = url
        self.port = 5984
        self.user = kwargs.get("username", None)
        self.pw = kwargs.get("password", None)
        self.url_string = "http://{}:{}".format(self.url, self.port)
        if log:
            self.log = log
        super(Couch, self).__init__(**kwargs)        
        if not self.con:
            raise ConnectionError("Connection failed for url {}".format(self.url_string))

    def connect(self, username=None, password=None, url="localhost", port=5984, **kw):
        if not username or not password or not url:
            self.log.warn("please supply username, password, and url")
            return None
        if not check_url(self.url_string):
            self.log.warn("No such url {}".format(self.url_string))
            return None
        self.con = couchdb.Server(url=self.url_string)
        self.log.info("Connected to server @{}".format(self.url_string))
        self.user = username
        self.pw = password

    def set_db(self, dbname):
        """Set database to use

        :param dbname: database name
        """
        try:
            self.db = self.con[dbname]
        except:
            return None

    def get_entry(self, name, field=None):
        """Retrieve entry from db for a given name, subset to field if
        that value is passed.

        :param name: unique name identifier (primary key, not the uuid)
        :param field: get 'field' of document, i.e. key in document dict
        """
        if not self._doc_type:
            return
        self.log.debug("retrieving field entry in field '{}' for name '{}'".format(field, name))
        if self.name_view.get(name, None) is None:
            self.log.warn("no field '{}' for name '{}'".format(field, name))
            return None
        doc = self._doc_type(**self.db.get(self.name_view.get(name)))
        if field:
            return doc[field]
        else:
            return doc

    def save(self, obj):
        """Save/update database object <obj>. If <obj> already exists
        and <update_fn> is defined, update will only take place if
        object has been modified

        :param obj: database object to save
        """
        if not self._update_fn:
            self.db.save(obj)
            self.log.info("Saving object {} with id {}".format(repr(obj), obj["_id"]))
        else:
            (new_obj, dbid) = self._update_fn(self.db, obj)
            if not new_obj is None:
                self.log.info("Saving object {} with id {}".format(repr(new_obj), new_obj["_id"]))
                self.db.save(new_obj)
            else:
                self.log.info("Object {} with id {} present and not in need of updating".format(repr(obj), dbid.id))


class GenoLogics(Database):
    def __init__(**kwargs):
        super(Couch, self).__init__(**kwargs)

    def connect(self, username=None, password=None):
        pass
