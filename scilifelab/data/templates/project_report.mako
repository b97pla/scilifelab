.. footer::

   \#\#\#Page\#\#\#

.. image:: ${sll_logo_small}
   :align: right

===================
Project status note
===================

${date}
---------------------------------

Description
^^^^^^^^^^^

This is a project status note containing an overview of the status of your project to date. If your You will receive a new status update

Project name
^^^^^^^^^^^^

${project_name} ${'({})'.format(customer_reference) if customer_reference not in ['', 'N/A'] else ''}

UPPNEX project id
^^^^^^^^^^^^^^^^^

${uppnex_project_id}

Sequence data directories
^^^^^^^^^^^^^^^^^^^^^^^^^

/proj/${uppnex_project_id}/INBOX/${project_name}

Samples
^^^^^^^

${sample_table}

Comments
^^^^^^^^

${finished}


Information
^^^^^^^^^^^

Naming conventions
""""""""""""""""""

The data is delivered in fastq format using Illumina 1.8 quality
scores. There will be one file for the forward reads and one file for
the reverse reads. More information on our naming conventions can be
found at
http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf.

Data access at UPPMAX
"""""""""""""""""""""

Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next
Generation sequence Cluster Storage, www.uppmax.uu.se), from which the
user can access it. If you have problems to access your data, please
contact SciLifeLab genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.
Information on how to access your data can be found at
http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf.

Acknowledgement
"""""""""""""""

In publications based on data from the work covered by this contract, the 
authors must acknowledge SciLifeLab, NGI and Uppmax: "The authors would 
like to acknowledge support from Science for Life Laboratory, the National 
Genomics Infrastructure, NGI, and Uppmax for providing assistance in massive 
parallel sequencing and computational infrastructure."
