import os
from mako.template import Template
from data import _files as data_files

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

def setup_empty_files():
    """Setup empty files"""
    for f in _empty_files():
        if not os.path.exists(os.path.join(filedir, f)):
            print "Preparing test: touching file {}".format(f)
            if not os.path.exists(os.path.dirname(os.path.join(filedir, f))):
                os.makedirs(os.path.dirname(os.path.join(filedir, f)))
            with open(os.path.join(filedir, f), "w") as fh:
                fh.write("")

def _empty_files():
    ## pre-casava result files
    files = ['data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_10_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_10_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_12_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_12_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_1_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_1_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_2_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_2_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_3_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_3_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_4_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_4_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_8_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_8_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_unmatched_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_unmatched_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.align_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.dup_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.hs_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.insert_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_17_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_17_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_19_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_19_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_5_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_5_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_7_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_7_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_unmatched_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_unmatched_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/3_120829_AA001AAAXX_barcode/3_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/4_120829_AA001AAAXX_barcode/4_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/5_120829_AA001AAAXX_barcode/5_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/6_120829_AA001AAAXX_barcode/6_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/7_120829_AA001AAAXX_barcode/7_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/8_120829_AA001AAAXX_barcode/8_120829_AA001AAAXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8.bam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8.sam',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8_1.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8_2.sai',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/1_120829_BB001BBBXX_barcode/1_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/2_120829_BB001BBBXX_barcode/2_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/3_120829_BB001BBBXX_barcode/3_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/4_120829_BB001BBBXX_barcode/4_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/5_120829_BB001BBBXX_barcode/5_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/6_120829_BB001BBBXX_barcode/6_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/7_120829_BB001BBBXX_barcode/7_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/8_120829_BB001BBBXX_barcode/8_120829_BB001BBBXX_nophix.bc_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix.filter_metrics',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/production/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/archive/120829_SN0001_0002_BB001BBBXX/run_info.yaml']

    ## Casava output 
    files.extend([
            'data/log/.placeholder',
            'data/projects/j_doe_00_01/README',
            'data/projects/j_doe_00_02/README',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort-dup.align_metrics',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort-dup.bam',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort-dup.dup_metrics',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort-dup.hs_metrics',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort-dup.insert_metrics',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5-sort.bam',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5_1_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_5_2_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort-dup.align_metrics',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort-dup.bam',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort-dup.dup_metrics',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort-dup.hs_metrics',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort-dup.insert_metrics',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2-sort.bam',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2_1_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/3_120914_BB002ABCXX_nophix_2_2_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort-dup.align_metrics',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort-dup.bam',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort-dup.dup_metrics',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort-dup.hs_metrics',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort-dup.insert_metrics',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9-sort.bam',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9_1_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/4_120914_BB002ABCXX_nophix_9_2_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort-dup.align_metrics',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort-dup.bam',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort-dup.dup_metrics',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort-dup.hs_metrics',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort-dup.insert_metrics',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1-sort.bam',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1_1_fastq.txt',
            'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/5_120914_BB002ABCXX_nophix_1_2_fastq.txt',
            'data/projects/j_doe_00_03/j_doe_00_03_git/config/post_process.yaml',
            ])
    files.extend(_add_casava_results(data_files()['data/archive/120924_SN0002_0003_CC003CCCXX/C003CCCXX.csv']))
    files.extend(_add_project_analyses())
    return files

## Add casava result files
def _add_casava_results(runinfo):
        ## Generate the sample files for casava 
    bcids = [7,2,5,7,3,8,4,1]
    i=0
    tmp = []
    for r in runinfo.split("\n"):
        v = r.split(",")
        if v[0] == "FCID":
            continue
        k = {'lane':v[1], 'name':v[2], 'sample_prj':v[9].replace("__", "."), 'sequence':v[4], 'barcode_id':bcids[i]}
        i = i + 1
        tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${name}_${sequence}_L00${lane}_R1_001.fastq").render(**k))
        tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${name}_${sequence}_L00${lane}_R2_001.fastq").render(**k))

        ## Add meta files to root folder
        file_types = ["-bcbb-command.txt","-bcbb-command.txtre","-bcbb.log","-post_process.yaml"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${name}${ext}").render(**k))

        ## Add meta files to root folder
        file_types = ["01_analysis_start.txt","02_process_lane.txt","03_remove_contaminants.txt","04_process_alignment.txt","05_merge_sample.txt","06_mark_duplicates_sample.txt","07_screen_sample_contaminants.txt","08_recalibrate_sample.txt","09_realign_sample.txt","10_variantcall.txt","11_detect_sv.txt","12_process_sample.txt","13_generate_bigwig.txt","14_write_metrics.txt", "bcbb_software_versions.txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${ext}").render(**k))
        tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${lane}_120924_CC003CCCXX.bc_metrics").render(**k))

        ## Add sort files to root folder
        file_types  = ["-dup-insert.pdf", "-dup-summary.aux", "-dup-summary.log", "-dup-summary.pdf", "-dup-summary.tex", "-dup.align_metrics", "-dup.bam", "-dup.dup_metrics", "-dup.hs_metrics", "-dup.insert_metrics", ".bam"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/${lane}_120924_CC003CCCXX_${barcode_id}_nophix-sort${ext}").render(**k))
            
        ## Add files to alignments folder
        file_types = ["-sort.bam",".bam",".sam","_1.sai","_1_fastq-fastq.bam","_2.sai"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/alignments/${lane}_120924_CC003CCCXX_${barcode_id}_nophix${ext}").render(**k))

        ## Add files to fastq_screen
        file_types = [".png", ".txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/fastq_screen/${lane}_120924_CC003CCCXX_${barcode_id}_nophix_1_fastq_screen${ext}").render(**k))

        ## Add files to fastqc
        file_types = ["fastqc_data.txt", "fastqc_report.html", "summary.txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/fastqc/${lane}_120924_CC003CCCXX_${barcode_id}_nophix-sort-dup_fastqc/${ext}").render(**k))
        file_types = ["_1_fastq.txt.gz", "_2_fastq.txt.gz", ".filter_metrics"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/production/${sample_prj}/${name}/120924_CC003CCCXX/nophix/${lane}_120924_CC003CCCXX_${barcode_id}_nophix${ext}").render(**k))
    return tmp

## Add project analyses
def _add_project_analyses():
    k={}
    ## final results data
    file_types = ["-dup-gatkrecal-realign-insert.pdf","-dup-gatkrecal-realign-summary.aux","-dup-gatkrecal-realign-summary.log","-dup-gatkrecal-realign-summary.pdf","-dup-gatkrecal-realign-summary.tex",
                  "-dup-gatkrecal-realign-variants-combined-phased-annotated.vcf","-dup-gatkrecal-realign-variants-combined-phased-effects.tsv","-dup-gatkrecal-realign-variants-combined-phased-effects.vcf",
                  "-dup-gatkrecal-realign-variants-combined-phased-effects.vcf.idx","-dup-gatkrecal-realign-variants-combined-phased.eval","-dup-gatkrecal-realign-variants-combined-phased.vcf",
                  "-dup-gatkrecal-realign-variants-combined-phased.vcf.eval_metrics","-dup-gatkrecal-realign-variants-combined-phased.vcf.idx","-dup-gatkrecal-realign-variants-combined.vcf",
                  "-dup-gatkrecal-realign-variants-combined.vcf.idx","-dup-gatkrecal-realign-variants-indel-filterINDEL.vcf","-dup-gatkrecal-realign-variants-indel-filterINDEL.vcf.idx",
                  "-dup-gatkrecal-realign-variants-indel.vcf","-dup-gatkrecal-realign-variants-indel.vcf.idx","-dup-gatkrecal-realign-variants-snp-SNPfilter.vcf",
                  "-dup-gatkrecal-realign-variants-snp-SNPfilter.vcf.idx","-dup-gatkrecal-realign-variants-snp.recal","-dup-gatkrecal-realign-variants-snp.tranches",
                  "-dup-gatkrecal-realign-variants-snp.vcf","-dup-gatkrecal-realign-variants-snp.vcf.idx","-dup-gatkrecal-realign-variants.vcf",
                  "-dup-gatkrecal-realign-variants.vcf.idx","-dup-gatkrecal-realign.align_metrics","-dup-gatkrecal-realign.bam",
                  "-dup-gatkrecal-realign.bam.bai","-dup-gatkrecal-realign.bigwig","-dup-gatkrecal-realign.hs_metrics","-dup-gatkrecal-realign.insert_metrics","-dup-gatkrecal-realign.pileup.gz",
                  "-dup-gatkrecal.bam","-dup-gatkrecal.bam.bai","-dup.bam","-dup.dup_metrics","-dup.recal",".bam"]
    tmp = []
    for x in file_types:
        k.update(ext=x)
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_7-sort${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_2-sort${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_7-sort${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_2-sort${ext}").render(**k))


    ## realign-split data
    file_types = ["-chr1-realign-subsetchr1.bam","-chr1-realign-subsetchr1.bam.bai","-chr1-realign.bam","-chr1-realign.intervals",
                  "-chr10-realign-subsetchr10.bam","-chr10-realign-subsetchr10.bam.bai","-chr10-realign.bam","-chr10-realign.intervals"]

    for x in file_types:
        k.update(ext=x)
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-split/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-split/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-split/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-split/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))

    file_types = ["-chr1-realign.bai","-chr10-realign.bai"]

    for x in file_types:
        k.update(ext=x)
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        
    ## variants-split data
    file_types = ["-chr1-variants.vcf","-chr1-variants.vcf.idx", "-chr10-variants.vcf", "-chr10-variants.vcf.idx"]
    for x in file_types:
        k.update(ext=x)
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-variants-split/1_120924_CC003CCCXX-sort-dup-gatkrecal-realign${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-variants-split/1_120924_CC003CCCXX-sort-dup-gatkrecal-realign${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-variants-split/1_120924_CC003CCCXX-sort-dup-gatkrecal-realign${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-variants-split/1_120924_CC003CCCXX-sort-dup-gatkrecal-realign${ext}").render(**k))

    file_types = ["-chr1-variants.vcf.idx","-chr10-variants.vcf.idx"]
    for x in file_types:
        k.update(ext=x)
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-variants-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_1/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-variants-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_2-sort-dup-gatkrecal-realign-variants-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))
        tmp.append(Template("data/projects/j_doe_00_04/intermediate/analysis_2/1_120924_CC003CCCXX_7-sort-dup-gatkrecal-realign-variants-split/tx/1_120924_CC003CCCXX-sort-dup-gatkrecal${ext}").render(**k))

    return tmp
