.. image:: ${sll_logo_small}
   :align: right

===================
Project status note
===================

SciLifeLab Stockholm
--------------------

${date}
-----------------

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
the reverse reads (if the run was a paired-end run). 

The naming of the files follow the convention:

  [LANE]_[DATE]_[FLOWCELL]_[SCILIFE NAME]_[READ].fastq.gz

Data access at UPPMAX
"""""""""""""""""""""

Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next
Generation sequence Cluster Storage, www.uppmax.uu.se), from which the
user can access it. You can find the data in the INBOX folder of the UPPNEX project, which
was created for you when your order was placed, e.g. 

  /proj/b2013000/INBOX/J.Doe_13_01

If you have problems to access your data, please
contact SciLifeLab genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.
 
Acknowledgement
"""""""""""""""

In publications based on data from the work covered by this contract, the 
authors must acknowledge SciLifeLab, NGI and Uppmax: "The authors would 
like to acknowledge support from Science for Life Laboratory, the National 
Genomics Infrastructure, NGI, and Uppmax for providing assistance in massive 
parallel sequencing and computational infrastructure."
