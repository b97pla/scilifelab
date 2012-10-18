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

    def __repr__(self):
        return "{}".format(self.__class__)

## From http://stackoverflow.com/questions/8780168/how-to-begin-writing-a-python-wrapper-around-another-wrapper
class Couch(Database):
    def __init__(self, log=None, **kwargs):
        self.db = None
        self.url = kwargs.get("url", None)
        self.port = 5984
        self.user = kwargs.get("username", None)
        self.pw = kwargs.get("password", None)
        self.url_string = "http://{}:{}".format(self.url, self.port)
        if log:
            self.log = log
        super(Couch, self).__init__(**kwargs)
        if not self.con:
            raise ConnectionError("Connection failed for url {}".format(self.url_string))
        
    def connect(self, username=None, password=None, url=None, port=5984):
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

        
class GenoLogics(Database):
    def __init__(**kwargs):
        super(Couch, self).__init__(**kwargs)

    def connect(self, username=None, password=None):
        pass
