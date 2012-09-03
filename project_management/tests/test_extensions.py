"""
Test extensions
"""

import os
from cement.core import handler
from pmtools.core.analysis import AnalysisController
from test_default import PmTest

class PmHsMetricsTest(PmTest):
    def test_1_hsmetrics(self):
        """Run hs metrics"""
        self.app = self.make_app(argv=['analysis', 'hs-metrics', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01', '-r', 'regionfile', '-n'], extensions=['pmtools.ext.ext_hs_metrics'])
        handler.register(AnalysisController)
        print self.app.argv
        self._run_app()

    def test_2_hsmetrics_empty(self):
        """Run hs metrics when no files present"""
        pass

    def test_3_hsmetrics_drmaa(self):
        """Run hs metrics over drmaa"""
        pass
