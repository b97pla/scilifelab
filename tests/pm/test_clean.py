"""
Test cleaning operations
"""

import os
import sys
import glob
from cement.core import handler
from cement.utils import shell
from test_default import PmTest, safe_makedir
from scilifelab.pm.core.project import ProjectController

flowcell = "120829_AA001AAAXX"
filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
intermediate = os.path.join(filedir, "data", "projects", "j_doe_00_02", "intermediate")
data = os.path.join(filedir, "data", "projects", "j_doe_00_02", "data")

class CleanTest(PmTest):
    INPUT_FILES = {'P1_106F_index6' : [
            '1_120829_AA001AAAXX_nophix_10_1_fastq.txt',
            '1_120829_AA001AAAXX_nophix_10_2_fastq.txt'],
                   'P1_107_index7' : [
            '1_120829_AA001AAAXX_nophix_12_1_fastq.txt.gz',
            '1_120829_AA001AAAXX_nophix_12_2_fastq.txt.gz']
                   }
    RESULT_FILES = {'P1_106F_index6': [
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign.bam.bai',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal.bai',
            '1_120829_AA001AAAXX_nophix_10-sort.bigwig',
            '1_120829_AA001AAAXX_nophix_10-sort-dup-target.crisp_pileup.gz',
            '1_120829_AA001AAAXX_nophix_10-sort-dup-target.pileup.gz',
            '1_120829_AA001AAAXX_nophix_10-sort-dup-target.crisp_pileup',
            '1_120829_AA001AAAXX_nophix_10-sort-dup-target.pileup',
            'alignments/1_120829_AA001AAAXX_nophix_10_1_fastq-fastq.bam',
            'alignments/1_120829_AA001AAAXX_nophix_10_2_fastq-fastq.bam',
            'alignments/1_120829_AA001AAAXX_nophix_10_1.sai',        
            'alignments/1_120829_AA001AAAXX_nophix_10_2.sai',        
            'alignments/1_120829_AA001AAAXX_nophix_10.sam',        
            'alignments/1_120829_AA001AAAXX_nophix_10.bam',        
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr1-realign-subsetchr1.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr1-realign-subsetchr1.bam.bai',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr1-realign.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr1-realign.intervals',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr2-realign-subsetchr2.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr2-realign-subsetchr2.bam.bai',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr2-realign.bam',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-chr2-realign.intervals',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-chr1-variants.vcf',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-chr2-variants.vcf',
            '1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_10-sort-gatkrecal-realign-chr3-variants.vcf'
            ],
                    'P1_107_index7' : [
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign.bam.bai',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal.bai',
            '1_120829_AA001AAAXX_nophix_12-sort.bigwig',
            '1_120829_AA001AAAXX_nophix_12-sort-dup-target.crisp_pileup.gz',
            '1_120829_AA001AAAXX_nophix_12-sort-dup-target.pileup.gz',
            '1_120829_AA001AAAXX_nophix_12-sort-dup-target.crisp_pileup',
            '1_120829_AA001AAAXX_nophix_12-sort-dup-target.pileup',
            'alignments/1_120829_AA001AAAXX_nophix_12_1_fastq-fastq.bam',
            'alignments/1_120829_AA001AAAXX_nophix_12_2_fastq-fastq.bam',
            'alignments/1_120829_AA001AAAXX_nophix_12_1.sai',        
            'alignments/1_120829_AA001AAAXX_nophix_12_2.sai',        
            'alignments/1_120829_AA001AAAXX_nophix_12.sam',        
            'alignments/1_120829_AA001AAAXX_nophix_12.bam',        
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr1-realign-subsetchr1.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr1-realign-subsetchr1.bam.bai',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr1-realign.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr1-realign.intervals',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr2-realign-subsetchr2.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr2-realign-subsetchr2.bam.bai',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr2-realign.bam',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-chr2-realign.intervals',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-chr1-variants.vcf',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-chr2-variants.vcf',
            '1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-variants-split/1_120829_AA001AAAXX_nophix_12-sort-gatkrecal-realign-chr3-variants.vcf'
            ]
                    }

    def setUp(self):
        super(CleanTest, self).setUp()
        ## Setup pre-casava results
        for k in self.RESULT_FILES.keys():
            for f in self.RESULT_FILES[k]:
                outfile = os.path.join(intermediate, flowcell, f)
                if not os.path.exists(os.path.dirname(outfile)):
                    safe_makedir(os.path.dirname(outfile))
                exit_code = shell.exec_cmd2(['touch', outfile])
        for k in self.INPUT_FILES.keys():
            for f in self.INPUT_FILES[k]:
                outfile = os.path.join(data, flowcell, '1_120829_AA001AAAXX_barcode', f)
                if not os.path.exists(os.path.dirname(outfile)):
                    safe_makedir(os.path.dirname(outfile))
                exit_code = shell.exec_cmd2(['touch', outfile])

        ## Setup casava results
        for k in self.RESULT_FILES.keys():
            for f in self.RESULT_FILES[k]:
                outfile = os.path.join(data, k, flowcell, f)
                if not os.path.exists(os.path.dirname(outfile)):
                    safe_makedir(os.path.dirname(outfile))
                exit_code = shell.exec_cmd2(['touch', outfile])
        for k in self.INPUT_FILES.keys():
            for f in self.INPUT_FILES[k]:
                outfile = os.path.join(data, k, flowcell, '1_120829_AA001AAAXX_barcode', f)
                if not os.path.exists(os.path.dirname(outfile)):
                    safe_makedir(os.path.dirname(outfile))
                exit_code = shell.exec_cmd2(['touch', outfile])

    def test_1_clean_dry(self):
        self.app = self.make_app(argv = ['project', 'clean', 'j_doe_00_02', '--pileup', '-n', '--intermediate', '--force'])
        handler.register(ProjectController)
        self._run_app()

    def test_1_clean(self):
        before = glob.glob(os.path.join(intermediate, "120829_AA001AAAXX", "*"))
        self.app = self.make_app(argv = ['project', 'clean', 'j_doe_00_02', '--pileup', '--intermediate', '--force'])
        handler.register(ProjectController)
        self._run_app()
        after = glob.glob(os.path.join(intermediate, "120829_AA001AAAXX", "*"))
        diff = [os.path.basename(x) for x in list(set(before).difference(set(after)))]
        self.eq(set(diff), set(['1_120829_AA001AAAXX_nophix_12-sort-dup-target.pileup.gz', '1_120829_AA001AAAXX_nophix_12-sort-dup-target.pileup', '1_120829_AA001AAAXX_nophix_12-sort-dup-target.crisp_pileup.gz', '1_120829_AA001AAAXX_nophix_12-sort-dup-target.crisp_pileup', '1_120829_AA001AAAXX_nophix_10-sort-dup-target.crisp_pileup.gz', '1_120829_AA001AAAXX_nophix_10-sort-dup-target.crisp_pileup', '1_120829_AA001AAAXX_nophix_10-sort-dup-target.pileup', '1_120829_AA001AAAXX_nophix_10-sort-dup-target.pileup.gz']))

    def test_2_clean_fastqbam(self):
        before = glob.glob(os.path.join(data, "P1_106F_index6/120829_AA001AAAXX/alignments", "*"))
        self.app = self.make_app(argv = ['project', 'clean', 'j_doe_00_02', '--data', '--fastqbam', '--force'])
        handler.register(ProjectController)
        self._run_app()
        after = glob.glob(os.path.join(data, "P1_106F_index6/120829_AA001AAAXX/alignments", "*"))
        diff = len(before) - len(after)
        self.eq(diff, 2)
        diff = [os.path.basename(x) for x in list(set(before).difference(set(after)))]
        self.eq(set(diff), set(['1_120829_AA001AAAXX_nophix_10_1_fastq-fastq.bam', '1_120829_AA001AAAXX_nophix_10_2_fastq-fastq.bam']))
