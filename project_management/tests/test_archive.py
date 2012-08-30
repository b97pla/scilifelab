"""
Test analysis subcontroller
"""

from cement.core import handler
from cement.utils import test
from test_default import PmTestApp, clean
from pmtools.controller.archive import ArchiveController

class PmArchiveTest(test.CementTestCase):
    app_class = PmTestApp

    def setUp(self):
        clean()
        
    def test_1_default(self):
        app = self.make_app(argv=['archive'])
        handler.register(ArchiveController)
        app.setup()
        app.run()
        app.close()

    def test_2_ls(self):
        app = self.make_app(argv=['archive', 'ls'])
        handler.register(ArchiveController)
        app.setup()
        app.run()
        app.close()
