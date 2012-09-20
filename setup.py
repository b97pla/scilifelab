#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

setup(name = "scilifelab",
      version = "0.2",
      author = "SciLife",
      author_email = "genomics@scilifelab.se",
      description = "Useful scripts for use at SciLife",
      license = "MIT",
      namespace_packages=["scilife"],
      ## Comment out scripts that aren't needed. Not all scripts have
      ## right permissions so either should be omitted or chmoded.
      scripts = ['scripts/runsizes.py',
                 'scripts/clean_up.py',
                 'scripts/datasets.py'
                 'scripts/deliver_454_data.sh',
                 'scripts/deliver_miseq.py',
                 'scripts/fastq_screen.py',
                 'scripts/fastq_unique.py',
                 'scripts/generate_run_info.py',
                 'scripts/miseq_data.py',
                 'scripts/pm',
                 'scripts/pm-init',
                 'scripts/project_data_delivery.py',
                 'scripts/project_management.py',
                 'scripts/quota_log.py',
                 'scripts/run_cmd.py',
                 'scripts/run_genome_filter.py',
                 'scripts/run_hs_metrics.py',
                 'scripts/runinfotools.py',
                 'scripts/runsizes.py',
                 'scripts/split_demultiplexed.py',
                 'scripts/storage_load.py',
                 'scripts/test_drmaa.py',
                 'scripts/uppmax_stats.py',
                 'scripts/bcbb_helpers/run_bcbb_pipeline.py',
                 'scripts/bcbb_helpers/report_to_gdocs.py',
                 'scripts/bcbb_helpers/process_run_info.py',
                 'scripts/bcbb_helpers/fastq_screen_wrapper.py',
                 'scripts/bcbb_helpers/software_versions.py',
                 ## For batch jobs: comment out scripts that aren't
                 ## needed - move to subfolder in scripts?
                 'batch/compress_fastq.sh',
                 'batch/compress_pbzip2.sh',
                 'batch/compress_pigz.sh',
                 'batch/illumina_run_slurm.sh',
                 'batch/md5sum.sh',
                 'batch/nosetests_run.sh',
                 'batch/pmd5deep.sh',
                 'batch/run_cmd.sh',
                 'batch/run_cmd_core.sh',
                 'batch/run_fastq_screen.sh',
                 'batch/run_fastqc.sh',
                 'batch/sbatch_job_tmpl.sh',
                 'batch/swestore.py',
                 ## Utilities for data delivery - move to subfolder in
                 ## scripts?
                 'utilities/assisted_delivery.py',
                 'utilities/bin_reads_by_quality.py',
                 'utilities/casava_data_delivery.py',
                 'utilities/count_barcodes.py',
                 'utilities/fastq_utils.py',
                 'utilities/fc_delivery_reports.py',
                 'utilities/generate_sbatch.py',
                 'utilities/project_status_note.py',
                 'utilities/read_illumina_summary_xml.py',
                 'utilities/runQC_to_statusdb.py',
                 'utilities/sample_delivery_note.py',
                 ],
      install_requires = [
          "bcbio-nextgen >= 0.2",
          "drmaa >= 0.5",
	  "sphinx >= 1.1.3",
	  "couchdb >= 0.8",
          "reportlab >= 2.5",
          "cement >= 2.0.2",
          "mock"
      ],
      test_suite = 'nose.collector',
      packages=['scilife'],
      ## package_data: install data/templates needed by modules
      package_data = {'scilife':['pm/templates/tpl/make/*', 'data/*']}
      )
os.system("git rev-parse --short --verify HEAD > ~/.scilifelab_version")
