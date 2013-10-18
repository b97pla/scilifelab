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
You will also receive a project status note summarizing the overall progress of your project.

- The 'Project summary' section contains a brief overview of the project and the mapping between the sample names submitted together with the samples and the internal SciLifeLab IDs assigned to the samples by us (Table 1).  
- The 'Sequencing summary' section contains a brief overview of the sequencing setup as well as more detailed information regarding the sequencing yields (Table 2a and Table 2b) and quality of the reads (Table 3)
- The 'Data delivery' sections contains basic information about the file formats, naming convention and paths to the delivered files.

Please don't hesitate to contact genomics_support@scilifelab.se if you have any questions or comments.

.. raw:: pdf

  PageBreak
  
  SetPageCounter 1

Project summary
^^^^^^^^^^^^^^^

:Project name: ${project_name}
:Customer reference: ${'{}'.format(customer_reference) if customer_reference not in ['', 'N/A'] else 'N/A'}
:UPPNEX project id: ${uppnex_project_id}
:Application: ${application}
##:Sequencing lanes ordered: ${lanes_ordered}
##:Number of samples: ${no_samples}

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

The yield for each lane on the flowcell is reported in Table 2a and the yield for each sample and lane is reported in Table 2b.

% if application == 'Finished library':
Since the sequenced library was prepared and pooled by the customer, we can only guarantee the total yield for the lane but not for individual samples in the pool. 
If the total yield (including undemultiplexed reads) does not reach ${lane_yield_cutoff}M in the first sequencing run, we will usually do a second run.
This may not be the case if, for example, the low yield is due to a problem with the library.
% endif
 
.. table:: **Table 2a.** The total number of read${' pair' if is_paired else ''}s (in millions) in each lane.${' An asterisk (*) next to a value indicates a yield less than what we guarantee.' if application == 'Finished library' else ''}

  .. class:: sample-table  

${lane_yield_table}

% if not application == 'Finished library':
Note that a sample may have been sequenced on multiple flowcells but this report only contains the results from one flowcell.
If the total yield for a sample does not reach the ordered amount in the first sequencing run, we will usually do a second run. 
This may not be the case if, for example, the sample failed the incoming QC checks. 
% endif 
 
.. table:: **Table 2a.** The number of read${' pair' if is_paired else ''}s (in millions) for each sample and lane

  .. class:: sample-table  

${sample_yield_table}

Quality
~~~~~~~

To asses the quality of the sequencing, we look at the sample-independent error rate (PhiX) and the sample-specific quality values reported by the instrument (Q30).
If any of these values fail to meet our requirements, we will investigate the problem and may re-run the affected sample(s). 

% if application == 'Finished library':
Note that when the libraries were prepared by the customer, we can only give guarantees for the sample-independent error rate (PhiX).
%endif

We guarantee a lane-wise PhiX error rate of ${phix_cutoff}% or below${' (averaged across the forward and reverse reads)' if is_paired else ''}.
 
The cutoff for the sample-specific error rate will depend on the sequencing platform, run mode and number of cycles. 
For ${run_setup} cycles on a Illumina ${instrument_version}${' in {} mode'.format('high output' if not run_mode == 'RapidRun' else 'rapid') if not instrument_version == 'MiSeq' else ''},
we require the Q30 to be ${sample_q30_cutoff}% or above.  

The Q30 value represents the percentage of bases with a quality value of 30 or above.
 
.. table:: **Table 3.** The reported sample-specific average quality, Q30 and lane-wise PhiX error rate. An asterisk (*) next to a value indicates a quality below our requirements.

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

Data from the sequencing will be uploaded to the UPPNEX (UPPMAX Next Generation sequence Cluster Storage, www.uppmax.uu.se), from which the user can access it. 

Sequence reads in FASTQ format, compressed with gzip, have been delivered to the INBOX of the UPPNEX account associated with your project. Data for each sample is placed under its own subfolder in the project path.
Under each sample folder, the data is further organized according to flowcell. Since this is a ${'paired' if is_paired else 'single'} end run, there will be ${'two files' if is_paired else 'one file'} 
for each sample and lane${', one for the forward read and one for the reverse read' if is_paired else ''}. The quality scale is Sanger / phred33 / Illumina 1.8+.

The FASTQ file names are constructed according to the convention: 

- [Lane]_[Start date]_${'[Flowcell position]' if not instrument_version == 'MiSeq' else ''}[Flowcell id]_[SciLifeLab ID]_[Read no].fastq.gz

Please refer to the `Sequencing FAQ`_ for detailed instructions on how to log in and access your files. If you have problems to access your data, please contact genomics_support@scilifelab.se. If you have
questions regarding UPPNEX, please contact support@uppmax.uu.se.

.. _Sequencing FAQ: http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf

.. raw:: pdf

  PageBreak

Acknowledgement
^^^^^^^^^^^^^^^

In publications based on data from the work covered by this contract, the 
authors must acknowledge SciLifeLab, NGI and Uppmax: "The authors would 
like to acknowledge support from Science for Life Laboratory, the National 
Genomics Infrastructure, NGI, and Uppmax for providing assistance in massive 
parallel sequencing and computational infrastructure."



