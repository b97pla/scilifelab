"""
Test analysis subcontroller
"""

import os
from cement.core import handler
from test_default import PmTest
from pmtools.core.analysis import AnalysisController

class PmAnalysisTest(PmTest):

    def test_1_ls(self):
        self.app = self.make_app(argv = ['analysis', 'ls'])
        handler.register(AnalysisController)
        self._run_app()
        self.eq(self.app._output_data['stdout'].getvalue(), '120829_SN0001_0001_AA001AAAXX\n120829_SN0001_0002_BB001BBBXX')

    def test_2_bcstats(self):
        self.app = self.make_app(argv = ['analysis','bcstats'])
        handler.register(AnalysisController)
        self._run_app()

    def test_3_deliver(self):
        """Test default delivery to project directory"""
        # self.app = self.make_app(argv = ['analysis', 'deliver', '120829_SN0001_0001_AA001AAAXX'])
        # handler.register(AnalysisController)
        # self._run_app()
        self.app = self.make_app(argv = ['analysis', 'deliver', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01'])
        handler.register(AnalysisController)
        self._run_app()

        #self.eq(" ".join(self.app._output_data['stdout'].getvalue().split()[0:12]), "1 Lane 1, J.Doe_00_01 A001AAAXX unknown Align_illumina SampleSheet 1 J.Doe_00_01 P1_101F_index1 ATCACG")
