"""
Test project subcontroller
"""
import os
from cement.core import handler
from test_default import PmTest
from pmtools.core.project import ProjectController
from pmtools.lib.runinfo import get_runinfo


class ProjectTest(PmTest):
    runinfo = os.path.join(os.path.curdir, "data", "archive", "120829_SN0001_0001_AA001AAAXX", "run_info.yaml")
    FASTQ_COMPRESS_FILES = [
        '1_120829_AA001AAAXX_nophix_10_1_fastq.txt',
        '1_120829_AA001AAAXX_nophix_10_2_fastq.txt',
        '1_120829_AA001AAAXX_nophix_12_1_fastq.txt',
        '1_120829_AA001AAAXX_nophix_12_2_fastq.txt',
        '1_120829_AA001AAAXX_nophix_1_1.fastq',
        '1_120829_AA001AAAXX_nophix_1_2.fastq',
        '1_120829_AA001AAAXX_nophix_2_1.fastq',
        '1_120829_AA001AAAXX_nophix_2_2.fastq',
        '1_120829_AA001AAAXX_nophix_3_1.fq',
        '1_120829_AA001AAAXX_nophix_3_2.fq',
        '1_120829_AA001AAAXX_nophix_4_1_fastq.txt',
        '1_120829_AA001AAAXX_nophix_4_2_fastq.txt',
        '1_120829_AA001AAAXX_nophix_8_1_fastq.txt',
        '1_120829_AA001AAAXX_nophix_8_2_fastq.txt',
        ]

    OUTPUT_FILES = []

    def setUp(self):
        super(ProjectTest, self).setUp()
        self.info = get_runinfo(self.runinfo, tab=True)
        print self.info
        ## Make 
            
    def test_1_project_deliver(self):
        self.app = self.make_app(argv = ['project', 'deliver'])
        handler.register(ProjectController)
        self._run_app()
        
    def test_2_project_data_delivery(self):
        pass

    def test_3_compress(self):
        """Test compression of project data"""
        pass

    def test_4_compress_distributed(self):
        """Test distributed compression of project data"""
        pass
