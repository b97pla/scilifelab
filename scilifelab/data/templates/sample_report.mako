.. image:: ${sll_logo_small}
   :align: right

======================
Raw data delivery note
======================

SciLifeLab Stockholm
-----------------------

${date}
---------------------------------

Project name
^^^^^^^^^^^^

${project_name} ${'({})'.format(customer_reference) if customer_reference not in ['', 'N/A'] else ''}

UPPNEX project id
^^^^^^^^^^^^^^^^^

${uppnex_project_id}

Flow cell id
^^^^^^^^^^^^

${FC_id}

Sequence data directory
^^^^^^^^^^^^^^^^^^^^^^^

/proj/${uppnex_project_id}/INBOX/${project_name}/${scilifelab_name}/${start_date}_${FC_id}

Sample
^^^^^^

${scilifelab_name} / ${customer_name}. Ordered amount: ${ordered_amount} million read${'{}'.format(' pair') if is_paired else ''}s.

Method
^^^^^^

Clustered on cBot and sequenced on ${instrument_version}${', {} run mode,'.format(run_mode) if run_mode not in ['', 'N/A'] else ''} 
according to manufacturer's instructions. Demultiplexing and conversion using
${casava_version}. The quality scale is Sanger / phred33 / Illumina 1.8+.

Results
^^^^^^^

${rounded_read_count} million read${'{}'.format(' pair') if is_paired else ''}s${' in lane with PhiX error rate {}%'.format(phix_error_rate) if phix_error_rate != 'N/A' else ''}. 
Average quality score ${avg_quality_score} (${pct_q30_bases}% bases >= Q30).

Comments
^^^^^^^^

${success}

