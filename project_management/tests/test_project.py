"""
Test project subcontroller
"""

import os
from cement.utils import test
from cement.core import backend
from pmtools import PmController, PmApp
from test_default import PmTestApp

class InitProjectTest(test.CementTestCase):
    app_class = PmTestApp
    
    def test_project_init(self):
        self.app.setup()
        self.ok(self.app.config.has_key('archive', 'root'))
        self.eq(self.app.config.get('analysis', 'root'), os.path.join(os.path.abspath(os.getcwd()), "data", "analysis"))


    def test_project_data_delivery(self):
        pass
