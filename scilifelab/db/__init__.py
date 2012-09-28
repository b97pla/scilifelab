import os
import sys
import couchdb

from scilifelab.utils.http import check_url

class Database(object):
    """Main database connection object for noSQL databases"""

    def __init__(self, **kwargs):
        self.con = None
        self.connect(**kwargs)

    def connect(self, **kwargs):
        pass

## From http://stackoverflow.com/questions/8780168/how-to-begin-writing-a-python-wrapper-around-another-wrapper
class Couch(Database):
    def __init__(self, **kwargs):
        self.db = None
        self.url = None
        super(Couch, self).__init__(**kwargs)

    def connect(self, username=None, password=None, url=None, port=5984):
        if not username or not password or not url:
            return None
        self.url = "http://{}:{}".format(url, port)
        if not check_url(self.url):
            return None
        self.con = couchdb.Server(url=self.url)

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
