import os
import unittest
from scilifelab.report.rl import make_example_note

filedir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))

class TestSampleDeliveryNote(unittest.TestCase):
    def setUp(self):
        pass

    def test_1_make_example_note(self):
        make_example_note(os.path.join(filedir, "test.pdf"))

