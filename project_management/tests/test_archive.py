"""
Test analysis subcontroller
"""
import os
import yaml
from cement.core import handler
from test_default import PmTest, PmTestOutputHandler
from pmtools.core.archive import ArchiveController
from pmtools.lib.runinfo import get_runinfo, dump_runinfo

class PmArchiveTest(PmTest):
    OUTPUT_FILES = []
    
    def test_1_default(self):
        """Test default archive function"""
        self.app = self.make_app(argv=['archive'])
        handler.register(ArchiveController)
        self._run_app()

    def test_2_ls(self):
        """Test archive list"""
        self.app = self.make_app(argv=['archive', 'ls'])
        handler.register(ArchiveController)
        self._run_app()
        self.eq(self.app._output_data['stdout'].getvalue(), '120829_SN0001_0001_AA001AAAXX\n120829_SN0001_0002_BB001BBBXX\n')
    
    def test_3_runinfo_default(self):
        """Test runinfo default list"""
        self.app = self.make_app(argv=['archive', 'runinfo'])
        handler.register(ArchiveController)
        self._run_app()
        self.eq(self.app._output_data['stderr'].getvalue(), 'Please provide flowcell id')

    def test_4_runinfo(self):
        """Test runinfo list"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX'])
        handler.register(ArchiveController)
        self._run_app()
        info = get_runinfo(os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml"))
        self.eq(self.app._output_data['stdout'].getvalue(), dump_runinfo(info))

    def test_5_runinfo_tab(self):
        """Test runinfo tab list"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-t'])
        handler.register(ArchiveController)
        self._run_app()
        runinfo_tab = get_runinfo(os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml"), tab=True)
        self.eq(self.app._output_data['stdout'].getvalue(), dump_runinfo(runinfo_tab, tab=True))

    def test_6_runinfo_list_projects(self):
        """Test runinfo list available projects"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-P'])
        handler.register(ArchiveController)
        self._run_app()
        self.eq(self.app._output_data['stdout'].getvalue(), 'available projects for flowcell 120829_SN0001_0001_AA001AAAXX:\n\tJ.Doe_00_01\n\tJ.Doe_00_02')
