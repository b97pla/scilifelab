"""
Test analysis subcontroller
"""

import os
import yaml
import shutil
from cement.core import handler
from cement.utils import shell
from test_default import PmTest
from scilifelab.pm.core.analysis import AnalysisController
from scilifelab.pm.utils.misc import walk

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
delivery_dir = os.path.abspath(os.path.join(filedir, "data", "projects", "j_doe_00_01", "data"))
intermediate_delivery_dir = os.path.abspath(os.path.join(filedir, "data", "projects", "j_doe_00_01", "intermediate"))

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

    def test_3_from_pre_to_casava_transfer(self):
        """Test casava transfer to project directory"""
        self.app = self.make_app(argv = ['analysis', 'transfer', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01', '--from_pre_casava'])
        handler.register(AnalysisController)
        self._run_app()
        res = shell.exec_cmd(["ls", "-1", os.path.join(delivery_dir, "P1_101F_index1", "120829_AA001AAAXX")])
        self.eq(['1_120829_AA001AAAXX_barcode', '1_120829_AA001AAAXX_nophix_1-sort-dup.align_metrics', '1_120829_AA001AAAXX_nophix_1-sort-dup.bam', '1_120829_AA001AAAXX_nophix_1-sort-dup.dup_metrics', '1_120829_AA001AAAXX_nophix_1-sort-dup.hs_metrics', '1_120829_AA001AAAXX_nophix_1-sort-dup.insert_metrics', '1_120829_AA001AAAXX_nophix_1-sort.bam', 'P1_101F_index1-bcbb-config.yaml', 'alignments'], res[0].split())
    
    def test_4_from_pre_to_pre_casava_transfer(self):
        """Test pre_casava transfer to project directory"""
        self.app = self.make_app(argv = ['analysis', 'transfer', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01', '--from_pre_casava', '--to_pre_casava', '--quiet'])
        handler.register(AnalysisController)
        self._run_app()
        ## Assert data output
        res = shell.exec_cmd(["ls", "-1", os.path.join(delivery_dir, "120829_AA001AAAXX", "1_120829_AA001AAAXX_barcode")])
        res_files = ['1_120829_AA001AAAXX_nophix_10_1_fastq.txt', '1_120829_AA001AAAXX_nophix_10_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_12_1_fastq.txt','1_120829_AA001AAAXX_nophix_12_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_1_1_fastq.txt','1_120829_AA001AAAXX_nophix_1_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_2_1_fastq.txt','1_120829_AA001AAAXX_nophix_2_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_3_1_fastq.txt','1_120829_AA001AAAXX_nophix_3_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_4_1_fastq.txt','1_120829_AA001AAAXX_nophix_4_2_fastq.txt',
                     '1_120829_AA001AAAXX_nophix_8_1_fastq.txt','1_120829_AA001AAAXX_nophix_8_2_fastq.txt']
        self.eq(set(res_files), set(res[0].split()))
        ## Assert intermediate delivery output 
        res = shell.exec_cmd(["ls", "-1", os.path.join(intermediate_delivery_dir, "120829_AA001AAAXX")])
        self.eq(['1_120829_AA001AAAXX_nophix_1-sort-dup.align_metrics','1_120829_AA001AAAXX_nophix_1-sort-dup.bam'], res[0].split()[0:2])
        self.eq(['1_120829_AA001AAAXX_nophix_8-sort-dup.insert_metrics','1_120829_AA001AAAXX_nophix_8-sort.bam', 'alignments'], res[0].split()[-3:])
        ## Assert pruned yaml file contents
        with open(os.path.join(delivery_dir, "120829_AA001AAAXX", "project_run_info.yaml")) as fh:
            runinfo_yaml = yaml.load(fh)
        self.eq(runinfo_yaml['details'][0]['multiplex'][0]['name'], 'P1_101F_index1')
        self.eq(runinfo_yaml['details'][0]['multiplex'][0]['description'], 'J.Doe_00_01_P1_101F_index1')
        self.eq(set(runinfo_yaml['details'][0]['multiplex'][0]['files']), set([os.path.join(delivery_dir,"120829_AA001AAAXX", "1_120829_AA001AAAXX_barcode", os.path.basename(x)) for x in ['1_120829_AA001AAAXX_nophix_1_1_fastq.txt','1_120829_AA001AAAXX_nophix_1_2_fastq.txt']]))
            
        

    def test_5_quiet(self):
        """Test pre_casava delivery to project directory with quiet flag"""
        self.app = self.make_app(argv = ['analysis', 'transfer', '120829_SN0001_0001_AA001AAAXX', '-p', 'J.Doe_00_01', '--from_pre_casava', '--to_pre_casava', '--quiet'])
        handler.register(AnalysisController)
        self._run_app()

