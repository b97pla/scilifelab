"""
Test analysis subcontroller
"""

import os
from cement.core import handler
from test_default import PmTest
from pmtools.controller.analysis import AnalysisController

class PmAnalysisTest(PmTest):

    def test_1_ls(self):
        self.app = self.make_app(argv = ['analysis', 'ls'])
        handler.register(AnalysisController)
        self._run_app()
        self.eq(self.app._output_data['stdout'], ['120829_SN0001_0001_AA001AAAXX\n120829_SN0001_0002_BB001BBBXX\n'])

    def test_2_bcstats(self):
        self.app = self.make_app(argv = ['analysis','bcstats'])
        handler.register(AnalysisController)
        self._run_app()

    def test_3_hsmetrics(self):
        self.app = self.make_app(argv = ['analysis', 'hs_metrics'])
        handler.register(AnalysisController)
        self._run_app()
        
