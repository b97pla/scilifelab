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

    def test_3_runinfo_default(self):
        app = self.make_app(argv=['archive', 'runinfo'])
        handler.register(ArchiveController)
        app.setup()
        app.run()
        app.close()

    def test_4_runinfo(self):
        app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX'])
        handler.register(ArchiveController)
        app.setup()
        app.run()
        app.close()

    def test_5_runinfo_tab(self):
        app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-t'])
        handler.register(ArchiveController)
        app.setup()
        app.run()
        app.close()
