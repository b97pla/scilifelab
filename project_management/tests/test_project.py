"""
Test initialization of projects
"""

# import os
# import unittest
from cement.utils import test
# from pm.cli.main import pm

# class TestApp(pm):
#     class Meta:
#         argv = []
#         config_files = []

#class InitProjectTest(unittest.TestCase):
class InitProjectTest(test.CementTestCase):
#    app_class = pm
    def setUp(self):
        print "setting up test"
        #       print pm

    def test_project_init(self):
        print "Testing project init"

    def test_project_data_delivery(self):
        print "Testing project data delivery"
