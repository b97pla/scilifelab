"""
Test analysis subcontroller
"""

from cement.core import handler
from cement.utils import test
from test_default import PmTestApp
from pmtools.controller.analysis import AnalysisController

class PmAnalysisTest(test.CementTestCase):
    app_class = PmTestApp

    def test_1_ls(self):
        app = self.make_app(argv = ['analysis', 'ls'])
        handler.register(AnalysisController)
        app.setup()
        app.run()
        app.close()

    def test_2_bcstats(self):
        app = self.make_app(argv = ['analysis','bcstats'])
        handler.register(AnalysisController)
        app.setup()
        app.run()
        app.close()
        
