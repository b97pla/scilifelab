"""
Test analysis subcontroller
"""

from cement.core import handler
from test_default import PmTest
from pmtools.controller.archive import ArchiveController

class PmArchiveTest(PmTest):
    OUTPUT_FILES = []
    
    def test_1_default(self):
        self.app = self.make_app(argv=['archive'])
        handler.register(ArchiveController)
        self._run_app()

    def test_2_ls(self):
        self.app = self.make_app(argv=['archive', 'ls'])
        handler.register(ArchiveController)
        self._run_app()

    def test_3_runinfo_default(self):
        self.app = self.make_app(argv=['archive', 'runinfo'])
        handler.register(ArchiveController)
        self._run_app()

    def test_4_runinfo(self):
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX'])
        handler.register(ArchiveController)
        self._run_app()

    def test_5_runinfo_tab(self):
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-t'])
        handler.register(ArchiveController)
        self._run_app()
