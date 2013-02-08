import os
import unittest
import couchdb
import socket
import logbook
import time

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

LOG = logbook.Logger(__name__)

class SciLifeTest(unittest.TestCase):
    """Test class for small unit tests that only test functions and do
    not rely on the existence of pipeline-generated data.
    """
    def setUp(self):
        pass

class SciLifeFullTest(unittest.TestCase):
    """Test class for tests that require that test data has been
    completely run through the pipeline. The class will setup test
    fastq files and run them through the pipeline once.
    """

    def setUp(self):
        pass


def has_couchdb_installation():
    # Try connecting to server
    has_couchdb = True
    try:
        server = couchdb.Server()
        dbstats = server.stats()
    except socket.error as e:
        has_couchdb = False
        LOG.info("To run db tests setup a local couchdb server at http://localhost:5984")
        time.sleep(1)
        pass
    return has_couchdb
