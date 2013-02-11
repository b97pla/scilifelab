import os
import unittest
import ConfigParser
import logbook
from scilifelab.db.statusdb import  _match_barcode_name_to_project_sample

from ..classes import has_couchdb_installation

filedir = os.path.abspath(__file__)
flowcells = ["120924_SN0002_0003_AC003CCCXX", "121015_SN0001_0002_BB002BBBXX"]
projects = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
project_dir = os.path.join(filedir, "data", "production")
has_couchdb = has_couchdb_installation()
DATABASES = ["samples-test", "projects-test", "flowcells-test"]

LOG = logbook.Logger(__name__)

class TestStatusdbFunctions(unittest.TestCase):
    """Tests for statusdb functions that don't require a couchdb connection"""
    def setUp(self):
        self.project_samples = {"P001_101_index1":{'customer_name':'1_index1'},
                                "P001_101":{'customer_name':'1_index1'},
                                "P001_102" : {'customer_name':'102_index2'},
                                "P001_103B" : {'customer_name':'103'},
                                "4_index4" : {'customer_name':'gnu4'},
                                "5":{'customer_name':'gnu5'},
                                "SAMPLE_6":{'customer_name':'sample_6'},
                                "SAMPLE_6A":{'customer_name':'sample_6a'},
                                "Gnu7_700bp" : {'customer_name':'Gnu7'},
                                "8_index" : {'customer_name':"8_gnu"},
                                }
                              
    def test_map_barcode_to_project_sample(self):
        """Test mapping barcode names to project sample names. Add
        failing cases here.
        """

        # Should fail since we don't want to automagically transform the hundreds
        bc = "001_1_index1"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True, True)
        self.assertEqual(None, res)

        # This should match P001_101 but fails for now
        bc = "1_index1"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True)
        self.assertEqual(None, res)

        bc = "P005_5B_index5"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True)
        self.assertEqual('5', res.get("sample_name"))

        bc = "P001_101_index1"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True)
        self.assertEqual('P001_101_index1', res.get("sample_name"))

        bc = "P001_102_index2"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True)
        self.assertEqual('P001_102', res.get("sample_name"))

        # This should match SAMPLE_6 but fails for now
        bc = "SAMPLE_6B_index6"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True)
        self.assertEqual(None, res)

        bc = "SAMPLE_6A_index6"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True, force=True)
        self.assertEqual('SAMPLE_6A', res.get("sample_name"))

        bc = "Gnu7_700bp"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True, force=True)
        self.assertEqual('Gnu7_700bp', res.get("sample_name"))

        # Will currently fail
        bc = "8"
        res = _match_barcode_name_to_project_sample(bc, self.project_samples, True, force=True)
        self.assertEqual(None, res)


