.. image:: ${sll_logo_small}
   :align: right

.. header:: Document version: 1.0.2

.. footer:: ###Page###

================================
Best practice analysis report
================================

:Project: ${project_name}
:Application: Sequence capture
:Date: ${date}

Project summary
---------------

The ${capturekit} target enrichment kit was used to prepare and
sequence target-enriched DNA. Sequence data was mapped to a reference,
and variants were called using a publically available software
pipeline.

The results are summarized in the following sections. Table 1 lists
the samples included in the analysis, whereas table 2 shows a brief
summary of sequencing yield and QC. Table 3 summarizes the sequence
capture performance. Finally, Table 4 summarizes variant calls and
their comparison to known variants in dbSNP.

Sample summary
^^^^^^^^^^^^^^

.. table:: **Table 1**. Listing of analyzed samples and corresponding barcode sequences.

${table_sample_summary}

Sequencing summary
^^^^^^^^^^^^^^^^^^

.. table:: **Table 2**. Sequencing summary. Columns show the following statistics: 1) Total: total number of reads, expressed with SI prefixes 2) Aligned: percentage of aligned reads 3) Pair duplicates: the percentage of pair duplicates in a sample 4) The inferred insert size 

${project_summary}

Target summary
^^^^^^^^^^^^^^

.. table:: **Table 3**. Target summary. Columns show the following statistics: 1) On target bases: The percentage of bases mapping on target 2) Mean target coverage: the mean coverage in the target region 3) X10x coverage targets: the percentage of targets with at least 10X coverage 4) Zero coverage targets: the percentage of targets with zero coverage

${project_target_summary}

dbSNP summary
^^^^^^^^^^^^^^

.. table:: **Table 4**. dbSNP summary.  Columns show the following statistics: 1) Total variations: the number of quality filtered called variants 2) In dbSNP: the fraction of called snps that are in dbSNP 3) Transition transversion all, dbSNP, novel: the transition transversion ratios in snps stratified by all, dbSNP and novel.


${project_dbsnp_summary}


Analysis pipeline
-----------------

Data has been analyzed using the bcbb package [1]_. After sequencing,
samples have been aligned to the reference genome and post-processed
for variant calling. Briefly, the steps include:

1. alignment with short read alignment program bwa [2]_
2. sorting and PCR duplicate removal with picard [3]_
3. calculation of mapping and enrichment statistics with picard and gatk [4]_
4. variant calling with gatk

See the following section for details on software and database
versions used.

Alignment
^^^^^^^^^

Alignments are stored in the Sequence Alignment/Map format [5]_ in
binary format in \*.bam files. The bam files have been indexed (\*.bai
files), allowing for rapid visualization with e.g. the Integrative
genomics viewer [6]_. The \*.bigwig files display coverages and can be
viewed as tracks in the IGV or a genome browser (UCSC genome browser
[7]_ or Ensembl [8]_). The \*metrics files contain picard summary
metrics. 


Variant calling
^^^^^^^^^^^^^^^

Variant calling follows a best practice procedure implemented at the
BROAD institute [9]_. Variant calling and evaluation is performed with
the GATK genotyping tools. Variant calling is performed with
UnifiedGenotyper, in which snps and indels are called. Called variants
are then evaluated with VariantEval, in which the called variants are
compared and assessed to common variants in dbsnp, hapmap, and
1000genomes. Variants are also filtered with respect to several
quality scores, such as variant confidence and strand bias. The end
result is a number of vcf files [10]_:

- \*-variants.vcf - raw variants
- \*-variants-snp.vcf - raw snp calls
- \*-variants-snp-SNPfilter.vcf - filtered snp calls
- \*-variants-indel.vcf - raw indel calls
- \*-variants-indel-filterINDEL.vcf - filtered indel calls


Filtering serves to filter out variants with low credibility and
quality. Additional filtering based on population-based allele
frequencies may therefore be required. The filtered calls are combined
and annotated with predicted damaging effects of the variants using
snpEff [11]_. The output from these operations is found in the
following files:

- \*-variants-combined-phased-annotated.vcf
- \*-variants-combined-phased-effects.vcf
- \*-variants-combined-phased-effects.tsv


Finally, project-summary.csv is a text file summary containing various
metrics, such as alignment statistics and data on the found
variations. For instructions on how to access data at UPPMAX, see our
FAQ [12]_.

Analysis settings
-----------------

Database versions
^^^^^^^^^^^^^^^^^

.. table:: **Table 6**. Database files used for variant calling. The 'train' databases are used by GATK to recalibrate and assign a well-calibrated probability to each variant call in a call set.

${database_versions_table}

Software versions
^^^^^^^^^^^^^^^^^

.. table:: **Table 7**. Software versions used by the pipeline.

${software_versions_table}


.. raw:: pdf

   PageBreak


References
----------

.. [1] URL: https://github.com/chapmanb/bcbb.

.. [2] Heng Li and Richard Durbin. "Fast and accurate short read
   alignment with Burrows-Wheeler transform". In: Bioinformatics 25.14
   (July 2009), pp. 1754-1760. URL: http://bioinformatics.oxford
   journals.org/cgi/content/abstract/25/14/1754.

.. [3] URL: http://picard.sourceforge.net/.

.. [4] URL: http://tinyurl.com/6oj5gj5

.. [5] Heng Li et al. "The Sequence Alignment/Map format and
   SAMtools". In: Bioinformatics (Oxford, England) 25.16 (Aug. 2009).
   PMID: 19505943, pp. 2078-2079. URL: http://www.ncbi.nlm.nih.gov/
   pubmed/19505943.

.. [6] James T Robinson et al. "Integrative genomics viewer". In: Nat
   Biotech 29.1 (Jan. 2011), pp. 24-26. URL:
   http://dx.doi.org/10.1038/nbt.1754.

.. [7] URL: genome.ucsc.edu.

.. [8] URL: www.ensembl.org.

.. [9] Mark A DePristo et al. "A framework for variation discovery and
   genotyping using next-generation DNA sequencing data". In: Nature
   Genetics 43.5 (May 2011). PMID: 21478889, pp. 491-498. URL:
   http://www.ncbi.nlm.nih.gov/pubmed/21478889.

.. [10] Petr Danacek et al. "The variant call format and VCFtools". In:
   Bioinformatics (Oxford, England) 27.15 (Aug. 2011). PMID: 21653522,
   pp. 2156-2158. URL: http://www.ncbi.nlm.nih.gov/pubmed/21653522.

.. [11] P Cingolani. snpEff: Variant effect prediction. 2012. URL:
   http://snpeff.sourceforge.net.

.. [12] URL: http://www.scilifelab.se/archive/pdf/tmp/SciLifeLab_Sequencing_FAQ.pdf

