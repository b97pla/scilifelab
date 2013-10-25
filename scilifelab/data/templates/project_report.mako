.. footer::

   \#\#\#Page\#\#\#

.. header::

  .. table::

${logo_table}

====================================================
${project_name} - Project Summary
====================================================

${date}
---------------------------------

Description
^^^^^^^^^^^

This report summarizes the progress of your project to date. If your samples are sequenced on multiple flowcells, you may receive
a project summary after each delivery. Note that this summary may contain information about sequencing results that you have not 
yet received. 

The report contains the following sections:

- `Project overview`_ contains a brief summary of the project parameters
- `Sequencing overview`_ contains an overview of the sequencing runs that have been carried out on the project samples
- `Sample overview`_ contains an overview of the progress of each sample in the project
- Sample details contains detailed information about the progress of each sample in the project

Acknowledgement
^^^^^^^^^^^^^^^

In publications based on data from the work covered by this contract, the 
authors must acknowledge SciLifeLab, NGI and Uppmax: "The authors would 
like to acknowledge support from Science for Life Laboratory, the National 
Genomics Infrastructure, NGI, and Uppmax for providing assistance in massive 
parallel sequencing and computational infrastructure."

Please don't hesitate to contact support@ngisweden.zendesk.com if you have any questions or comments.

.. raw:: pdf

  PageBreak
  
  SetPageCounter 1

Project overview
^^^^^^^^^^^^^^^^

:Project name: ${project_name}
:Customer reference: ${'{}'.format(customer_reference) if customer_reference not in ['', 'N/A'] else 'N/A'}
:Application: ${application}
${':Reads per sample ordered: {}'.format(m_ordered) if m_ordered else ''}
${':Sequencing lanes ordered: {}'.format(lanes_ordered) if lanes_ordered else ''}
${':Number of samples: {}'.format(no_samples) if no_samples else ''}
:UPPNEX project id: ${uppnex_project_id}
:Project path: /proj/${uppnex_project_id}/INBOX/${project_name}

Sequencing overview
^^^^^^^^^^^^^^^^^^^

Table 1 contains a summary of the sequencing lanes that have been run to date for the ${project_name} project.

.. table:: **Table 1** The sequencing lanes where samples from the project have been sequenced

  .. class:: sample-table  

${flowcell_table}

Sample overview
^^^^^^^^^^^^^^^

Table 2 contains an overview of the progress of each sample in the project. Only reads that have passed our quality criteria are included in the read counts below.

.. table:: **Table 2** Summary of the progress of each sample in the project.

  .. class:: sample-table  

${sample_table}

.. raw:: pdf

  PageBreak

Data delivery
^^^^^^^^^^^^^

Data from the sequencing are uploaded to the UPPNEX (UPPMAX Next
Generation sequence Cluster Storage, www.uppmax.uu.se), from which the
user can access it. If you have problems to access your data, please
contact SciLifeLab genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.
Detailed information on how to access your data can be found in the `Sequencing FAQ`_ 
as well as in the sample delivery note that you received together with each delivery.

.. _Sequencing FAQ: http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf
