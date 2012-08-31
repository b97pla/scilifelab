"""
Test analysis subcontroller
"""
import os
import yaml
from cement.core import handler
from test_default import PmTest, PmTestOutputHandler
from pmtools.controller.archive import ArchiveController
from pmtools.lib.runinfo import runinfo_to_tab

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
        self.eq(self.app._output_data['stdout'], ['120829_SN0001_0001_AA001AAAXX\n120829_SN0001_0002_BB001BBBXX\n'])
    
    def test_3_runinfo_default(self):
        """Test runinfo default list"""
        self.app = self.make_app(argv=['archive', 'runinfo'])
        handler.register(ArchiveController)
        self._run_app()
        self.eq(self.app._output_data['stderr'], ['Please provide flowcell id'])

    def test_4_runinfo(self):
        """Test runinfo list"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX'])
        handler.register(ArchiveController)
        self._run_app()
        with open(os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml")) as fh:
            runinfo_data = yaml.load(fh)
        self.eq(self.app._output_data['stdout'][0], runinfo_data)


    def test_5_runinfo_tab(self):
        """Test runinfo tab list"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-t'])
        handler.register(ArchiveController)
        self._run_app()
        with open(os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml")) as fh:
            runinfo_data = yaml.load(fh)
        runinfo_data_tab = runinfo_to_tab(runinfo_data)
        self.eq(self.app._output_data['stdout'][0], runinfo_data_tab)

    def test_6_runinfo_list_projects(self):
        """Test runinfo list available projects"""
        self.app = self.make_app(argv=['archive', 'runinfo', '120829_SN0001_0001_AA001AAAXX', '-P'])
        handler.register(ArchiveController)
        self._run_app()
        self.eq(self.app._output_data['stdout'], ['available projects for flowcell 120829_SN0001_0001_AA001AAAXX:\n\tJ.Doe_00_01\n\tJ.Doe_00_02'])
