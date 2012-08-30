"""
Test analysis subcontroller
"""

from cement.core import handler
from cement.utils import test
from pmtools import PmApp
from test_default import PmTestApp
from pmtools.controller.analysis import AnalysisController

class PmAnalysisTest(test.CementTestCase):
    app_class = PmTestApp

    def test_1_deliver_fastq(self):
        app = self.make_app(argv = ['analysis', '-h'])
        handler.register(AnalysisController)
        self.app.setup()
        self.app.run()
        self.app.close()
