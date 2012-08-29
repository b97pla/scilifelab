"""
Test project subcontroller
"""

import os
from cement.utils import test
from pmtools import PmController, PmApp
from test_default import PmTestApp

class InitProjectTest(test.CementTestCase):
    app_class = PmTestApp

    def setUp(self):
        super(InitProjectTest, self).setUp()
        self.app = PmTestApp()
    
    def test_1_project_init(self):
        self.app.setup()
        
    def test_2_project_data_delivery(self):
        pass
