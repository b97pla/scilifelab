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
^^^^^^^

.. table:: **Sample names** Name conversion between SciLifeLab sample ID and the submitted sample ID 

${sample_name_table}

Yield
^^^^^

.. table:: **Sample yield** The number of million reads or read pairs resulting from the sequencing  

${sample_yield_table}

Quality
^^^^^^^

.. table:: **Sample quality** The reported sample quality and lane-wise PhiX error rate 

${sample_quality_table}
    
Method
^^^^^^

Sequenced on Illumina ${instrument_version}${' in {} mode'.format('high output' if not run_mode == 'RapidRun' else 'rapid') if not instrument_version == 'MiSeq' else ''} according to manufacturer's instructions. 
The sequencing setup was ${run_setup}bp. 
Basecalling was performed on the instrument using ${basecall_software}${' v{}'.format(basecaller_version) if basecaller_version else ''}. 
Bcl to fastq conversion and demultiplexing was performed using ${demultiplex_software}${' v{}'.format(demultiplex_version) if demultiplex_version else ''} from the CASAVA software suite.
The quality scale is Sanger / phred33 / Illumina 1.8+.
