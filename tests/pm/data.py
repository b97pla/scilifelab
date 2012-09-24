from mako.template import Template

files = {
    ### Old-school project analysis data
    'data/analysis/120829_SN0001_0001_AA001AAAXX/1_120829_AA001AAAXX_barcode/1_120829_AA001AAAXX_nophix.bc_metrics' : 
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
    ### Casava result files after delivery
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
""",
    ### Casava flowcell information
    'data/archive/120924_SN0002_0003_CC003CCCXX/C003CCCXX.csv':
        """FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
C003CCCXX,1,P001_101_index3,hg19,TGACCA,J__Doe_00_04,N,R1,NN,J__Doe_00_04
C003CCCXX,1,P001_102_index6,hg19,ACAGTG,J__Doe_00_04,N,R1,NN,J__Doe_00_04
C003CCCXX,2,P002_101_index3,hg19,TGACCA,J__Doe_00_05,N,R1,NN,J__Doe_00_05
C003CCCXX,2,P002_102_index6,hg19,ACAGTG,J__Doe_00_05,N,R1,NN,J__Doe_00_05
C003CCCXX,2,P002_103_index8,hg19,TGGTCA,J__Doe_00_05,N,R1,NN,J__Doe_00_05
C003CCCXX,2,P003_101_index1,hg19,AGTGCG,J__Doe_00_06,N,R1,NN,J__Doe_00_06
C003CCCXX,2,P003_102_index2,hg19,TGTGCG,J__Doe_00_06,N,R1,NN,J__Doe_00_06
C003CCCXX,2,P003_103_index6,hg19,CGTTAA,J__Doe_00_06,N,R1,NN,J__Doe_00_06""",
    ### Casava analysis data structures
    'data/analysis/120924_SN0002_0003_CC003CCCXX/1_120924_CC003CCCXX.bc_metrics':
        """7       22463443        TGACCA  P001_101_index3
2       63340036        ACAGTG  P001_102_index6
unmatched       2326234 Undetermined    lane1
""",
    'data/analysis/120924_SN0002_0003_CC003CCCXX/2_120924_CC003CCCXX.bc_metrics':
        """5       2246343        TGACCA  P002_101_index3
7       6334036        ACAGTG  P002_102_index6
3       4495853        TGGTCA  P002_103_index8
8       479491        AGTGCG  P003_101_index1
4       9316653        TGTGCG  P003_102_index2
1       7108259        CGTTAA  P003_103_index6
unmatched       3946195 Undetermined    lane2
"""}

## Generate the sample files for casava 
bcids = [7,2,5,7,3,8,4,1]
for r in files['data/archive/120924_SN0002_0003_CC003CCCXX/C003CCCXX.csv'].split("\n"):
    v = r.split(",")
    if v[0] == "FCID":
        continue
    k = {'lane':v[1], 'name':v[2], 'sample_prj':v[8].replace("__", "."), 'sequence':v[4]}
    files["data/analysis/{}/{}/120924_CC003CCCXX/{}-bcbb-config.yaml".format(v[5].replace("__", "."), v[2], v[2])] = Template("""details:
- analysis: Standard
  description: Lane ${lane}, ${sample_prj}
  flowcell_id: CC003CCCXX
  genome_build: hg19
  lane: '${lane}'
  multiplex:
  - analysis: Align_standard_seqcap
    barcode_id: 5
    barcode_type: SampleSheet
    description: ${sample_prj}_${name}
    files:
    - ${name}_${sequence}_L00${lane}_R1_001.fastq
    - ${name}_${sequence}_L00${lane}_R2_001.fastq
    genome_build: hg19
    genomes_filter_out: phix
    name: ${name}
    sample_prj: ${sample_prj}
    sequence: ${sequence}
fc_date: '120924'
fc_name: CC003CCCXX
""").render(**k)

