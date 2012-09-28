import os
import unittest
import ConfigParser
from scilifelab.db.statusdb import SampleRunMetricsConnection

class TestDbConnection(unittest.TestCase):
    def setUp(self):
        if not os.path.exists(os.path.join(os.getenv("HOME"), "dbcon.ini")):
            self.url = None
        else:
            config = ConfigParser.ConfigParser()
            config.readfp(open(os.path.join(os.getenv("HOME"), "dbcon.ini")))
            print config["couchdb"]

    def test_1_connection():
        print "hello"
