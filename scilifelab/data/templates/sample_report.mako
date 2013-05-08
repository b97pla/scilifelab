.. image:: ${sll_logo_small}
   :align: right

======================
Raw data delivery note
======================

SciLifeLab Stockholm
-----------------------

${date}
---------------------------------

Description
^^^^^^^^^^^

This is a raw data delivery note containing detailed information about the sequencing 
of your samples on one flowcell. If your samples have been sequenced on multiple flowcells, you will receive one delivery note for each flowcell.
You will also receive a project sample note summarizing the progress of your project.

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

/proj/${uppnex_project_id}/INBOX/${project_name}/[SciLifeLab ID]/${start_date}_${FC_id}

Samples
^^^^^^

${sample_table}

Method
^^^^^^

Clustered on ${'cBot' if not clustered == "OnBoardClustering" else 'board'} and 
sequenced on ${instrument_version} in ${'high output' if not run_mode == 'RapidRun' else 'rapid'} mode 
according to manufacturer's instructions. The sequencing setup was ${run_setup}bp. 
Basecalling was performed on instrument with ${basecall_software} v${basecaller_version}. 
Demultiplexing and fastq conversion were done using ${demultiplex_version}.
The quality scale is Sanger / phred33 / Illumina 1.8+.
