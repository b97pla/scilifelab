.. footer::

   \#\#\#Page\#\#\#

.. table::

${logo_table}

======================
Raw data delivery note
======================

${date}
---------------------------------

Description
^^^^^^^^^^^

This is a raw data delivery note containing detailed information about the sequencing 
of samples on one flowcell. If your samples have been sequenced on multiple flowcells, you will receive separate delivery notes for each flowcell.
You will also receive a project status note summarizing the progress of your project.

- **Table 1** shows the translation between the internal SciLifeLab ID given to each sample and the sample ID submitted together with the samples.
- **Table 2** shows a breakdown of the sequencing yield, in million read${' pair' if is_paired else ''}s for each sample and lane. Note that if a sample has been 
  run on multiple lanes on this flowcell and/or has been multiplexed using different barcode sequences, the sample will have one entry for each lane 
  and/or barcode.
- **Table 3** shows the reported quality for each sample.

Please don't hesitate to contact genomics_support@scilifelab.se if you have any questions or comments.

.. raw:: pdf

  PageBreak
  
  SetPageCounter 1

Project summary
^^^^^^^^^^^^^^^

:Project name: ${project_name}
:Customer reference: ${'{}'.format(customer_reference) if customer_reference not in ['', 'N/A'] else 'N/A'}
:UPPNEX project id: ${uppnex_project_id}

Upon arrival, the samples are assigned internal SciLifeLab IDs which are used in our documentation and processing. Table 1 specifies the translation between 
sample names submitted by the customer and the internal SciLifeLab ID. 

.. table:: **Table 1.** Name conversion between SciLifeLab sample ID and the submitted sample ID

  .. class:: sample-table  

${sample_name_table}

.. raw:: pdf

  PageBreak

Sequencing summary
^^^^^^^^^^^^^^^^^^

:Instrument type: Illumina ${instrument_version}
:Instrument ID: ${instrument_id}
:Sequencing setup: ${run_setup} cycles
:Run mode: ${'{}'.format('High Output' if not run_mode == 'RapidRun' else 'Rapid') if not instrument_version == 'MiSeq' else 'N/A'}
:Flowcell ID: ${FC_id}
:Flowcell position: ${FC_position}
:Sequencing start: ${start_date}


Method
~~~~~~

Sequenced on Illumina ${instrument_version}${' in {} mode'.format('high output' if not run_mode == 'RapidRun' else 'rapid') if not instrument_version == 'MiSeq' else ''} according to manufacturer's instructions. 
The sequencing setup was ${run_setup}bp. 
Basecalling was performed on the instrument using ${basecall_software}${' v{}'.format(basecaller_version) if basecaller_version else ''}. Bcl to fastq conversion and demultiplexing was performed using ${demultiplex_software}${' v{}'.format(demultiplex_version) if demultiplex_version else ''} from the CASAVA software suite.


Yield
~~~~~

The sequencing yield is reported for each sample and lane and/or barcode. 
 
.. table:: **Table 2.** The number of million read${' pair' if is_paired else ''}s resulting from this sequencing run

  .. class:: sample-table  

${sample_yield_table}


Quality
~~~~~~~

The quality reported in Table 3 is the percentage of bases with a quality score of 30 or above (Q>=30), the average quality score over all bases (Avg Q) and the
lane-wise sequencing error rate, estimated from sequencing and mapping of a PhiX spike-in (PhiX error rate).   

.. table:: **Table 3.** The reported sample quality and lane-wise PhiX error rate

  .. class:: sample-table 

${sample_quality_table}

.. raw:: pdf

  PageBreak

Data delivery
^^^^^^^^^^^^^

:Project path: /proj/${uppnex_project_id}/INBOX/${project_name}
:Sample path: [SciLifeLab ID]/${start_date}_${FC_position}${FC_id}
:File format: FASTQ
:Quality scale: Sanger / phred33 / Illumina 1.8+

Sequence reads in FASTQ format have been delivered to the INBOX of the UPPNEX account associated with your project. Data for each sample is placed under its own subfolder in the project path.
Under each sample folder, the data is further organized according to flowcell. Since this is a ${'paired' if is_paired else 'single'} end run, there will be ${'two files' if is_paired else 'one file'} 
for each sample and lane${', one for the forward read and one for the reverse read' if is_paired else ''}. The quality scale is Sanger / phred33 / Illumina 1.8+.

The FASTQ file names are constructed according to the convention: 

  [Lane]_[Start date]_${'[Flowcell position]' if not instrument_version == 'MiSeq' else ''}[Flowcell id]_[SciLifeLab ID]_[Read no].fastq.gz

Please refer to the `Sequencing FAQ`_ for detailed instructions on how to log in and access your files.

.. _Sequencing FAQ: http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf


