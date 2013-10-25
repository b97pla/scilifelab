
Sample details - ${scilifeid}
------------------------------------

Sample summary
^^^^^^^^^^^^^^

:SciLifeLab ID: ${scilifeid}
:Submitted ID: ${customerid}
${":Sample type: {}".format(sample_type) if sample_type else ""}
${":Ordered amount: {} M".format(m_ordered) if m_ordered else ""}
:Sequenced amount: ${"{} M".format(m_reads_sequenced) if m_reads_sequenced else "N/A"}
:Status: ${"Pass" if status else "N/A"}

Arrival QC
^^^^^^^^^^

${':Sample received: {}'.format(incoming_qc_start_date.strftime('%Y-%m-%d')) if incoming_qc_start_date else ''}
${':Arrival QC finished: {}'.format(incoming_qc_finish_date.strftime('%Y-%m-%d')) if incoming_qc_finish_date else ''}
:Arrival QC status: ${'N/A' if incoming_qc_status is None else '{}'.format('Pass' if incoming_qc_status else 'Fail')}
 
Library prep
^^^^^^^^^^^^

.. table::

  .. class:: sample-table  

${prep_table}

Sequencing runs
^^^^^^^^^^^^^^^

.. table::

  .. class:: sample-table  

${sequencing_table}
