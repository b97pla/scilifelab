"""
Test analysis subcontroller
"""

import os
from cement.core import handler
from test_default import PmTest
from pmtools.controller.analysis import AnalysisController

class PmAnalysisTest(PmTest):
    OUTPUT_FILES = []

    def test_1_ls(self):
        self.app = self.make_app(argv = ['analysis', 'ls'])
        handler.register(AnalysisController)
        self._run_app()

    def test_2_bcstats(self):
        self.app = self.make_app(argv = ['analysis','bcstats'])
        handler.register(AnalysisController)
        self._run_app()
        
