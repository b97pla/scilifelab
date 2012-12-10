================================
Analysis report sequence capture
================================
:Info: See <http://docutils.sf.net/rst.html> for introductory docs.

:abstract: The Agilent SureSelect All Exon 50mb target enrichment kit
           was used to prepare and sequence DNA enriched for exon sequence.
	   Sequence data was mapped to a reference, and variants were called
	   using a publically available software pipeline.

Project summary
---------------

Table ${tableref} shows a brief project summary. In the case of many
samples, the table is split into several tables.

.. table:: Project summary. Samples are listed by column. Each row
.. shows a summary statistics: 1) Total: total number of reads 2)
.. Aligned: number of aligned reads 3) Pair duplicates: the percentage
.. of pair duplicates in a sample 4) The inferred insert size 5) On
.. target bases: The percentage of bases mapping on target 6) Mean
.. target coverage: the mean coverage in the target region 7) X10x
.. coverage targets: the percentage of targets with at least 10X
.. coverage 8) Zero coverage targets: the percentage of targets with
.. zero coverage 9) Total variations: the number of quality filtered
.. called variants 10) In dbSNP: the fraction of called snps that are
.. in dbSNP 11) Transition transversion all, dbSNP, novel: the
.. transition transversion ratios in snps stratified by all, dbSNP and
.. novel.

${projectsummarytable}


Sample summary
--------------

Table 2 summarizes the samples that have been included in this
analysis.

.. table:: Sample summary, showing on which flowcell and which lane a
.. given sample has been run. The barcode sequence is included for
.. reference. If several runs have been performed on the same samples,
.. the analysis pipeline merges the results. Therefore, the final
.. result file name prefixes are not identical to the delivered data.
.. See section 4.3 for more information about naming conventions.

${samplesummarytable}

Data access at UPPMAX
---------------------

Data is delivered to a folder in the INBOX of your UPPMAX project
account. The delivery folder will typically be named DATE_IDENTIFIER,
where the identifier relates to the flowcell or analysis that has been
performed. The data can be accessed via ssh, scp or sftp. The
following tutorial describes how to access your data (change PRJID to
your project name).

Connecting to UPPMAX
^^^^^^^^^^^^^^^^^^^^

First, open an ssh session and change directory to your project folder:

.. code:: shell

    ssh yourUserName@biologin.uppmax.uu.se
    cd /proj/PRJID


We strongly recommend that you move the data from your INBOX to the
private directory in order to prevent unauthorized access:

.. code:: shell

   mv INBOX/DATE_IDENTIFIER private/

Enter your data directory and list your files:

.. code:: shell

   cd private/DATE_IDENTIFIER/
   ls

To check the size of your files, try:

.. code:: shell

   ls -lh
   du -h

Copying data from UPPMAX
^^^^^^^^^^^^^^^^^^^^^^^^

You can copy some or all of your data from UPPMAX to your local
computer. BEWARE: the output generated from analysis of
next-generation sequencing data can easily count in the hundreds of
gigabytes. In particular the mapping files (file extension .bam) are
large, and you may want to limit your download to a subset of the
data. Due to the massive data amounts, also keep in mind that data
transfer can take quite some time.

For example, if you want to copy only the variant call files (*.vcf
files) to your local computer, issue the following command (written on
one single line):

.. code:: shell

   rsync -avh -progress
   yourUserName@biologin.uppmax.uu.se:/proj/PRJID/private/DATE_IDENTIFIER/*.vcf .

Naming conventions
^^^^^^^^^^^^^^^^^^

Table 2 lists the samples included in the analysis. Delivered raw data
files have the prefix LANE DATE - FLOWCELLID, followed by a sample
identifier. The downstream analysis may induce a renaming of files, in
particular if the analysis includes data from several flowcells. The
naming convention in this case is that each sample is identified by
the file name prefix LANE SAMPLEID ANALYSISDATE ANALYSISID. Here,
ANALYSISDATE corresponds to the date of analysis, and ANALYSISID is a
unique identifier for the analysis.

Analysis pipeline
-----------------

Data has been analyzed using the bcbb package [2]. After sequencing,
samples have been aligned to the reference genome and post-processed
for variant calling. Briefly, the steps include:

1. alignment with short read alignment program bwa [9]
2. sorting and PCR duplicate removal with picard [1]
3. calculation of mapping and enrichment statistics with picard and gatk [8]
4. variant calling with gatk

Alignment
^^^^^^^^^

Alignments are stored in the Sequence Alignment/Map format [10] in
binary format in \*.bam files. The bam files have been indexed (\*.bai
files), allowing for rapid visualization with e.g. the Integrative
genomics viewer [11]. The \*.bigwig files display coverages and can be
viewed as tracks in the IGV or a genome browser (UCSC genome
browser[3] or Ensembl[4]). The \*metrics files contain picard summary
metrics.


Variant calling
^^^^^^^^^^^^^^^

Variant calling follows a best practice procedure implemented at the
BROAD institute [5]. Variant calling and evaluation is performed with
the GATK genotyping tools. Variant calling is performed with
UnifiedGenotyper, in which snps and indels are called. Called variants
are then evaluated with VariantEval, in which the called variants are
compared and assessed to common variants in dbsnp, hapmap, and
1000genomes. Variants are also filtered with respect to several
quality scores, such as variant confidence and strand bias. The end
result is a number of vcf files [7], where the most important are:

- \*-variants.vcf - raw variants
- \*-variants-snp.vcf - raw snp calls
- \*-variants-snp-SNPfilter.vcf - filtered snp calls
- \*-variants-indel.vcf - raw indel calls
- \*-variants-indel-filterINDEL.vcf - filtered indel calls

In addition to variant calling and filtering, the damaging effects of
the variants are predicted using snpEff [6]. The results are
summarized in the files snpEff_genes.txt and snpEff_summary.html.
project-summary.csv is a text file summary of the found variations.

References
----------

[1] URL: http://picard.sourceforge.net/.
[2] URL: https://github.com/chapmanb/bcbb.
[3] URL: genome.ucsc.edu.
[4] URL: www.ensembl.org.
[5] URL: http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3.
[6] P Cingolani. snpEff: Variant effect prediction. 2012. URL: http://snpeff.sourceforge.net.
[7] Petr Danacek et al. "The variant call format and VCFtools". In: Bioinformatics (Oxford, England) 27.15 (Aug. 2011). PMID: 21653522, pp. 2156-2158. URL: http://www.ncbi.nlm.nih.gov/pubmed/21653522.
[8] Mark A DePristo et al. "A framework for variation discovery and genotyping using next-generation DNA sequencing data". In: Nature Genetics 43.5 (May 2011). PMID: 21478889, pp. 491-498. URL: http://www.ncbi.nlm.nih.gov/pubmed/21478889.
[9] Heng Li and Richard Durbin. "Fast and accurate short read alignment with Burrows-Wheeler transform".
In: Bioinformatics 25.14 (July 2009), pp. 1754-1760. URL: http://bioinformatics.oxford
journals.org/cgi/content/abstract/25/14/1754.
[10] Heng Li et al. "The Sequence Alignment/Map format and SAMtools". In: Bioinformatics (Oxford, England)
25.16 (Aug. 2009). PMID: 19505943, pp. 2078-2079. URL: http://www.ncbi.nlm.nih.gov/
pubmed/19505943.
[11] James T Robinson et al. "Integrative genomics viewer". In: Nat Biotech 29.1 (Jan. 2011), pp. 24-26. URL:
http://dx.doi.org/10.1038/nbt.1754.
