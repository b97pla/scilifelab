
Sample summary - ${sample.scilife_name} (${sample.customer_name})
-----------------------------------------------------------------

:Arrival QC: ${sample.initial_QC_status or "N/A"}
:Library QC: ${sample.library_prep_status()}
:Average fragment size (customer): ${sample.customer_average_fragment_length}
:M reads sequenced: ${sample.m_reads_sequenced}
:Status: ${sample.status or "N/A"}

${sample.library_prep_table()}

${sample.sample_run_table(project)}

${sample.sample_delivery_table()}

##/proj/${uppnex_project_id}/INBOX/${project_name}/${scilifelab_name}/${start_date}_${FC_id}


##Clustered using ${clustering_method} and sequenced on ${sequencing_platform} (${sequencing_software}) with a ${sequencing_setup} setup in ${sequencing_mode} mode.
##Bcl to Fastq conversion was performed using bcl2Fastq v1.8.3 from the CASAVA software suite. The quality scale is Sanger / phred33 / Illumina 1.8+.
