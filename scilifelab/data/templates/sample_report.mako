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

${scilifelab_name} / ${customer_name}. Ordered amount: ${ordered_amount} million paired reads.

Method
^^^^^^

Clustered on cBot and sequenced on HiSeq 2000 according to
manufacturer's instructions. Base conversion using OLB v1.9,
demultiplexed and converted to fastq using CASAVA v1.8. The quality
scale is Sanger / phred33 / Illumina 1.8+. 

Results
^^^^^^^

${rounded_read_count} million reads in lane with PhiX error rate
${phix_error_rate}. Average quality score ${avg_quality_score}.

Comments
^^^^^^^^

${success}

