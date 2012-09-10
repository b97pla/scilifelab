"""
Test analysis subcontroller
"""

import os
import shutil
from cement.core import handler
from test_default import PmTest
from pmtools.core.analysis import AnalysisController
from pmtools.utils.misc import walk

delivery_dir = os.path.abspath(os.path.join(os.path.curdir, "data", "projects", "j_doe_00_01", "data"))

class PmAnalysisTest(PmTest):

    def setUp(self):
        super(PmAnalysisTest, self).setUp()
        if os.path.exists(delivery_dir):
            flist = walk(os.path.join(delivery_dir))
            for x in flist:
                os.unlink(x)
            shutil.rmtree(delivery_dir)
        
    def test_1_ls(self):
        self.app = self.make_app(argv = ['analysis', 'ls'])
        handler.register(AnalysisController)
        self._run_app()
        self.eq(self.app._output_data['stdout'].getvalue(), '120829_SN0001_0001_AA001AAAXX\n120829_SN0001_0002_BB001BBBXX')

    def test_2_bcstats(self):
        self.app = self.make_app(argv = ['analysis','bcstats'])
        handler.register(AnalysisController)
        self._run_app()

    def test_3_casava_delivery(self):
        """Test casava delivery to project directory"""
        self.app = self.make_app(argv = ['analysis', 'deliver', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01'])
        handler.register(AnalysisController)
        self._run_app()

    def test_4_pre_casava_delivery(self):
        """Test pre_casava delivery to project directory"""
        self.app = self.make_app(argv = ['analysis', 'deliver', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01', '--pre_casava'])
        handler.register(AnalysisController)
        self._run_app()
        
