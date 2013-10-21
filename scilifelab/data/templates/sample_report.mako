
.. table::

${logo_table}

${date} - Sample status
---------------------------------

===============================================================
${scilifeid} ${" - {}".format(customerid) if customerid else ""}
===============================================================

Description
^^^^^^^^^^^

This is a sample status note containing detailed information about the progress of one sample. 
The aim of this report is to provide an overall picture of the progress of the sample and the 
information within reflects the status at the date at the top. Note that this report may contain
information about sequencing runs that you haven't yet received any data for. 
 
If your sample is sequenced on multiple flowcells, you may receive
an updated sample status note after each delivery.

Please don't hesitate to contact genomics_support@scilifelab.se if you have any questions or comments.

Arrival QC
^^^^^^^^^^

${':Sample received: {}'.format(incoming_qc_start_date.strftime('%Y-%m-%d')) if incoming_qc_start_date else ''}
${':Arrival QC finished: {}'.format(incoming_qc_finish_date.strftime('%Y-%m-%d')) if incoming_qc_finish_date else ''}
:Arrival QC status: ${'N/A' if incoming_qc_status is None else '{}'.format('Passed' if incoming_qc_status else 'Failed')}
 
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
