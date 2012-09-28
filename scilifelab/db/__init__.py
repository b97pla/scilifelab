import os
import sys
import couchdb

from cement.core import backend

LOG = backend.minimal_logger("db")

from scilifelab.utils.http import check_url

class Database(object):
    """Main database connection object for noSQL databases"""

    def __init__(self, **kwargs):
        self.con = None
        self.log = LOG
        self.connect(**kwargs)

    def connect(self, **kwargs):
        pass

## From http://stackoverflow.com/questions/8780168/how-to-begin-writing-a-python-wrapper-around-another-wrapper
class Couch(Database):
    def __init__(self, log=None, **kwargs):
        self.db = None
        self.url = None
        if log:
            self.log = log
        super(Couch, self).__init__(**kwargs)

    def connect(self, username=None, password=None, url=None, port=5984):
        if not username or not password or not url:
            self.log.warn("please supply username, password, and url")
            return None
        self.url = "http://{}:{}".format(url, port)
        if not check_url(self.url):
            self.log.warn("No such url {}".format(self.url))
            return None
        self.con = couchdb.Server(url=self.url)
        self.log.info("Connected to server @{}".format(self.url))

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
