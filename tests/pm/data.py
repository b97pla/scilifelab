files = {'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix.bc_metrics' : 
         """1       19756915
2       18724985
3       21948744
4       17394069
8      26933322
10      18014097
12	23252366
unmatched       9289601
""",
         'data/analysis/120829_SN0001_0001_AA001AAAXX/2_120829_AA001AAAXX_barcode/2_120829_AA001AAAXX_nophix.bc_metrics' :
             """5       19235231
7       10232523
17       2194
19       17125
unmatched       22001241
""",
         'data/archive/120829_SN0001_0001_AA001AAAXX/run_info.yaml':
"""- analysis: Align_illumina
  description: Lane 1, J.Doe_00_01
  flowcell_id: A001AAAXX
  genome_build: unknown
  lane: '1'
  multiplex:
  - barcode_id: 1
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_101F_index1
    sample_prj: J.Doe_00_01
    sequence: ATCACG
  - barcode_id: 2
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_102F_index2
    sample_prj: J.Doe_00_01
    sequence: CGATGT
  - barcode_id: 3
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_103_index3
    sample_prj: J.Doe_00_01
    sequence: TTAGGC
  - barcode_id: 4
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_104F_index4
    sample_prj: J.Doe_00_01
    sequence: TGACCA
  - barcode_id: 8
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_105F_index5
    sample_prj: J.Doe_00_01
    sequence: ACAGTG
  - barcode_id: 10
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_106F_index6
    sample_prj: J.Doe_00_01
    sequence: GCCAAT
  - barcode_id: 12
    barcode_type: SampleSheet
    genome_build: unknown
    name: P1_107_index7
    sample_prj: J.Doe_00_01
    sequence: CAGATC
- analysis: Align_illumina
  description: Lane 2, J.Doe_00_02
  flowcell_id: A001AAAXX
  genome_build: unknown
  lane: '2'
  multiplex:
  - barcode_id: 5
    barcode_type: SampleSheet
    genome_build: unknown
    name: P2_101_index19a
    sample_prj: J.Doe_00_02
    sequence: ATCACG
  - barcode_id: 7
    barcode_type: SampleSheet
    genome_build: unknown
    name: P2_102_index12a
    sample_prj: J.Doe_00_02
    sequence: CGATGT
  - barcode_id: 17
    barcode_type: SampleSheet
    genome_build: unknown
    name: P2_103_index3a
    sample_prj: J.Doe_00_02
    sequence: TTAGGC
  - barcode_id: 19
    barcode_type: SampleSheet
    genome_build: unknown
    name: P2_104_index4a
    sample_prj: J.Doe_00_02
    sequence: TGACCA
""",
         'data/projects/j_doe_00_03/data/P000_101/120914_BB002ABCXX/P000_101-bcbb-config.yaml':
"""details:
- analysis: Align_illumina
  description: Lane 3, J.Doe_00_03
  flowcell_id: B002ABCXX
  genome_build: unknown
  lane: '3'
  multiplex:
  - analysis: Align_illumina
    barcode_id: 5
    barcode_type: SampleSheet
    description: J.Doe_00_03_P000_101
    files: [./3_120914_BB002ABCXX_nophix_5_2_fastq.txt,
            ./3_120914_BB002ABCXX_nophix_5_1_fastq.txt]
    name: P000_101
    sample_prj: J.Doe_00_03
    sequence: ATCACG
""",
         'data/projects/j_doe_00_03/data/P000_102/120914_BB002ABCXX/P000_102-bcbb-config.yaml':
"""details:
- analysis: Align_illumina
  description: Lane 3, J.Doe_00_03
  flowcell_id: B002ABCXX
  genome_build: unknown
  lane: '3'
  multiplex:
  - analysis: Align_illumina
    barcode_id: 2
    barcode_type: SampleSheet
    description: J.Doe_00_03_P000_102
    files: [./3_120914_BB002ABCXX_nophix_2_2_fastq.txt,
            ./3_120914_BB002ABCXX_nophix_2_1_fastq.txt]
    name: P000_102
    sample_prj: J.Doe_00_03
    sequence: ATCACG
""",
         'data/projects/j_doe_00_03/data/P000_103/120914_BB002ABCXX/P000_103-bcbb-config.yaml':
"""details:
- analysis: Align_illumina
  description: Lane 4, J.Doe_00_03
  flowcell_id: B002ABCXX
  genome_build: unknown
  lane: '4'
  multiplex:
  - analysis: Align_illumina
    barcode_id: 9
    barcode_type: SampleSheet
    description: J.Doe_00_03_P000_103
    files: [./4_120914_BB002ABCXX_nophix_9_2_fastq.txt,
            ./4_120914_BB002ABCXX_nophix_9_1_fastq.txt]
    name: P000_103
    sample_prj: J.Doe_00_03
    sequence: ATCACG
""",
         'data/projects/j_doe_00_03/data/P000_104F/120914_BB002ABCXX/P000_104F-bcbb-config.yaml':
"""details:
- analysis: Align_illumina
  description: Lane 5, J.Doe_00_03
  flowcell_id: B002ABCXX
  genome_build: unknown
  lane: '5'
  multiplex:
  - analysis: Align_illumina
    barcode_id: 1
    barcode_type: SampleSheet
    description: J.Doe_00_01_P1_101F_index1
    files: [./5_120914_BB002ABCXX_nophix_1_2_fastq.txt,
            ./5_120914_BB002ABCXX_nophix_1_1_fastq.txt]
    name: P000_104F_index1
    sample_prj: J.Doe_00_03
    sequence: ATCACG
"""
         }
