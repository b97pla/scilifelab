galaxy_config: universe_wsgi.ini
program:
  bowtie: bowtie
  samtools: samtools
  bwa: bwa
  maq: maq
  ucsc_bigwig: wigToBigWig
  picard: $PICARD_HOME
  gatk: $GATK_HOME
  fastqc: fastqc
  pdflatex: pdflatex
  ps2pdf: ps2pdf
  barcode: barcode_sort_trim.py
  snpEff: $SNPEFF_HOME
algorithm:
  aligner: bwa
  demultiplexed: true
  max_errors: 2
  num_cores: 8
  stringency: high
  quality_format: Illumina
  platform: illumina
  recalibrate: true
  snpcall: true 
  bc_mismatch: 1
  bc_allow_indels: false
  bc_read: 1
  bc_position: 3
  bc_illumina_no_trailing: true
  bc_offset: 1
  java_memory: 3g
  save_diskspace: true
  screen_contaminants: false
  filter_phix: false
  upload_fastq: true

distributed:
  run_process_program_locally: true
  cluster_platform: slurm
  num_workers: 16
  cores_per_host: 1
  rabbitmq_vhost: bionextgen
  platform_args: -p core -A a2010002 -t 2:00:00


analysis:
  store_dir: ${store_dir}
  base_dir: ${base_dir}
  towig_script: echo
  distributed_process_program: run_bcbb_pipeline.py
  process_program: automated_initial_analysis.py

# configuration algorithm changes for custom post-processing steps
custom_algorithms:
  'Broad SNP':
    ref_ext: '-broad'
    aligner: maq
    recalibrate: true
    snpcall: true
    dbsnp: variation/dbsnp_132.vcf
    hybrid_bait: bait_list
    hybrid_target: target_list
  'SNP calling':
    aligner: bwa
    recalibrate: true
    snpcall: true
    dbsnp: variation/dbsnp_132.vcf
  'Minimal':
    aligner: ""
  'Align_illumina':
    aligner: bwa
    recalibrate: false
    snpcall: false
    java_memory: 3g
    quality_format: Illumina
    upload_fastq: true
    hybrid_target: 
    hybrid_bait: 
  'Align_standard':
    aligner: bwa
    recalibrate: false
    snpcall: false
    java_memory: 3g
    quality_format: Standard
  'Align_standard_seqcap':
    aligner: bwa
    recalibrate: true
    snpcall: true
    java_memory: 3g
    quality_format: Standard
    upload_fastq: true
    dbsnp: ${dbsnp}
    train_hapmap: ${hapmap}
    train_1000g_omni: ${omni}
    train_indels: ${mills}

#    hybrid_target: 
#    hybrid_bait: 

resources:
  ucsc_bigwig:
    memory: 3g
  gatk:
    cores: 8
  bwa:
    cores: 8
  mosaik:
    cores: 1
  bowtie:
    cores: 8
  tophat:
    cores: 1
  cufflinks:
    cores: 1
