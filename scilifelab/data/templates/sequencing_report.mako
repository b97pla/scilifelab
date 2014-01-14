
Project details - ${project.get("project_name")}
-------------------------------------------------


:Application: ${application}
${':Reads per sample ordered: {}'.format(m_ordered) if m_ordered else ''}
${':Sequencing lanes ordered: {}'.format(lanes_ordered) if lanes_ordered else ''}
${':Number of samples: {}'.format(no_samples) if no_samples else ''}
:UPPNEX project id: ${uppnex_project_id}
:Project path: /proj/${uppnex_project_id}/INBOX/${project_name}

Upon arrival, the samples are assigned internal SciLifeLab IDs which are used in our documentation and processing. Table 1 specifies the translation between 
sample names submitted by the customer and the internal SciLifeLab ID. 

.. table:: **Table 1.** Name conversion between SciLifeLab sample ID and the submitted sample ID

  .. class:: sample-table  

${sample_name_table}

% if type == 'applications':
Note that your project is being handled by the Genomics Applications team, which typically means that some 
aspects of the project contains features that are under development. This may have consequences on e.g. the
expected time required to finish the project or the quality of the data that we guarantee. Please don't hesitate to
contact support@ngisweden.zendesk.com if you have any questions.
% endif 

Sequencing summary
~~~~~~~~~~~~~~~~~~

Table 2 contains a summary of the sequencing lanes that have been run to date on project samples.

.. table:: **Table 2** The sequencing lanes where samples from the project have been sequenced

  .. class:: sample-table  

${flowcell_table}

Sample overview
~~~~~~~~~~~~~~~

Table 3 contains an overview of the progress of each sample in the project. Only reads that have passed our quality criteria are included in the read counts below.
The status column indicates the status of each sample. Possible values are 'Finished', 'In progress' and 'Aborted'. 

.. table:: **Table 3** Summary of the progress of each sample in the project.

  .. class:: sample-table  

${sample_table}

.. raw:: pdf

  PageBreak

