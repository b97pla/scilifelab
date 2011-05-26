!SLIDE 
# Blue Collar Bioinformatics NGS pipeline #

!SLIDE bullets incremental center
# Original Idea #
![chapman](chapman.jpg)
### Brad Chapman ###
#### Massachussets General Hospital ####
#### **Very** active bioinformatics developer ####
#### twitter:@chapmanb github:chapmanb ####

!SLIDE bullets incremental center
# Contributors #
### Roman Valls Guimera ###
#### Science for Life Laboratory Stockholm ####
#### Not to confuse with Sci-Fi ####
#### http://scilifelab.se ####
#### twitter:@braincode github:brainstorm ####

!SLIDE center
![github.png](github.png)

!SLIDE center bullets incremental
# Vision #
* Illumina sequencer ends run
* Convert to standard FastQ format
* Align against reference genome
* SNP calling, variant eval & effect
* Summary PDF report

## "Automate, fix, improve, automate" ##

!SLIDE bullets incremental transition-fade
# Objectives #
* Provide an automatic **preliminar** analysis
* Transform data, ready to use on UPPMAX: SAM, BAM, BAM-sorted, bigwig...
* Perform adhoc analysis with Galaxy
* Perform further adhoc analysis outside Galaxy

!SLIDE small center
# Automation #
    "(...) It [processing] does take a bit, but that's why it's
    nice to have it automated, it'll all go
    off in the middle of night
    without having to think about it."

                            -- Brad Chapman


!SLIDE bullets incremental transition=scrollUp
## Good coders code, great **reuse** ##
### NGS analysis ###
 * Picard
 * GATK
 * bowtie
 * bwa
 * tophat
 * samtools

!SLIDE bullets incremental transition=scrollUp
## NGS analysis (cont) ###
 * snpEff
 * fastx toolkit
 * FastQC 
 * matrix2png

!SLIDE bullets incremental transition=scrollUp
## Python modules ##
 * Biopython
 * pysam
 * mako
 * PyYAML
 * amqplib
 * logbook
 * nosetests

!SLIDE bullets incremental transition=scrollUp
## Processing infrastructure ##
 * RabbitMQ
 * LaTeX

!SLIDE bullets incremental transition=scrollUp
## Optional software for generating report graphs ##
 * R with ggplot2, plyr, sqldf libraries.
 * rpy2



!SLIDE center small
# Collaborative software development #
![softdevel.png](softdevel.png)

!SLIDE center small
# Flexibility vs Complexity #

![groucho.jpg](groucho.jpg)

##"Those are my principles, if you don't like them, I've others"##
#### -- Groucho Marx ####

!SLIDE center small 
# "Talk is cheap, show me the code" #
#### -- Linus Torvalds (Linux Kernel hacker) ####

!SLIDE center small commandline incremental


    $ cd bcbb/nextgen/scripts && find
    ./picard_gatk_recalibrate.py
    ./picard_maq_recalibrate.py
    ./analyze_quality_recal.py
    ./automated_initial_analysis.py
    ./bam_to_wiggle.py
    ./analyze_finished_sqn.py
    ./utils/resort_bam_karyotype.py
    ./utils/broad_redo_analysis.py
    ./utils/convert_samplesheet_config.py
    ./utils/sort_gatk_intervals.py
    ./solexa_qseq_to_fastq.py
    ./upload_to_galaxy.py
    ./store_finished_sqn.py
    ./picard_sam_to_bam.py
    ./variant_effects.py
    ./align_summary_report.py
    ./gatk_variant_eval.py
    ./monthly_billing_report.py
    ./barcode_sort_trim.py
    ./illumina_finished_msg.py
    ./gatk_genotyper.py

    $ find . -iname "*.py" | xargs wc -l 
    (...)
    4264 totalt

