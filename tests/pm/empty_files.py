from mako.template import Template
from data import files as data_files

def files():
    ## pre-casava result files
    files = ['data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_10_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_10_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_12_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_12_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_1_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_1_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_2_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_2_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_3_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_3_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_4_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_4_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_8_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_8_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_unmatched_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix_unmatched_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_1-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_10-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_12-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_2-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_3-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_4-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.align_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.dup_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.hs_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort-dup.insert_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_nophix_8-sort.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_17_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_17_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_19_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_19_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_5_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_5_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_7_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_7_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_unmatched_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix_unmatched_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/3_120829_AA001AAAXX_barcode/3_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/4_120829_AA001AAAXX_barcode/4_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/5_120829_AA001AAAXX_barcode/5_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/6_120829_AA001AAAXX_barcode/6_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/7_120829_AA001AAAXX_barcode/7_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/8_120829_AA001AAAXX_barcode/8_120829_AA001AAAXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_10_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_12_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_1_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_2_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_3_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_4_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8.bam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8.sam',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8_1.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/alignments/1_120829_AA001AAAXX_nophix_8_2.sai',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/1_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/2_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/3_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/4_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/5_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/6_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/7_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0001_AA001AAAXX/nophix/8_120829_AA001AAAXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/1_120829_BB001BBBXX_barcode/1_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/2_120829_BB001BBBXX_barcode/2_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/3_120829_BB001BBBXX_barcode/3_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/4_120829_BB001BBBXX_barcode/4_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/5_120829_BB001BBBXX_barcode/5_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/6_120829_BB001BBBXX_barcode/6_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/7_120829_BB001BBBXX_barcode/7_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/8_120829_BB001BBBXX_barcode/8_120829_BB001BBBXX_nophix.bc_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/1_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/2_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/3_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/4_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/5_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/6_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/7_120829_BB001BBXX_nophix_2_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix.filter_metrics',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix_1_fastq.txt',
             'data/analysis/120829_SN0001_0002_BB001BBBXX/nophix/8_120829_BB001BBXX_nophix_2_fastq.txt',
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
    files.extend(add_casava_results(data_files()['data/archive/120924_SN0002_0003_CC003CCCXX/C003CCCXX.csv']))
    return files

## Add casava result files
def add_casava_results(runinfo):
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
        tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${name}_${sequence}_L00${lane}_R1_001.fastq").render(**k))
        tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${name}_${sequence}_L00${lane}_R2_001.fastq").render(**k))
        file_types = ["-bcbb-command.txt","-bcbb-command.txtre","-bcbb.log","-post_process.yaml"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${name}${ext}").render(**k))
        file_types = ["01_analysis_start.txt","02_process_lane.txt","03_remove_contaminants.txt","04_process_alignment.txt","05_merge_sample.txt","06_mark_duplicates_sample.txt","07_screen_sample_contaminants.txt","08_recalibrate_sample.txt","09_realign_sample.txt","10_variantcall.txt","11_detect_sv.txt","12_process_sample.txt","13_generate_bigwig.txt","14_write_metrics.txt", "bcbb_software_versions.txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${ext}").render(**k))
        tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${lane}_120924_CC003CCCXX.bc_metrics").render(**k))
        file_types  = ["-insert.pdf", "-summary.aux", "-summary.log", "-summary.pdf", "-summary.tex", ".align_metrics", ".bam", ".dup_metrics", ".hs_metrics", ".insert_metrics"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/${lane}_120924_CC003CCCXX_${barcode_id}_nophix-sort-dup${ext}").render(**k))
        file_types = ["-sort.bam",".bam",".sam","_1.sai","_1_fastq-fastq.bam","_2.sai"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/alignments/${lane}_120924_CC003CCCXX_${barcode_id}_nophix${ext}").render(**k))
        file_types = [".png", ".txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/fastq_screen/${lane}_120924_CC003CCCXX_${barcode_id}_nophix_1_fastq_screen${ext}").render(**k))
        file_types = ["fastqc_data.txt", "fastqc_report.html", "summary.txt"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/fastqc/${lane}_120924_CC003CCCXX_${barcode_id}_nophix-sort-dup_fastqc/${ext}").render(**k))
        file_types = ["_1_fastq.txt.gz", "_2_fastq.txt.gz", ".filter_metrics"]
        for x in file_types:
            k.update(ext=x)
            tmp.append(Template("data/analysis/${sample_prj}/${name}/120924_CC003CCCXX/nophix/${lane}_120924_CC003CCCXX_${barcode_id}_nophix${ext}").render(**k))
    return tmp





