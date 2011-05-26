!SLIDE
# How do we barcode ?

!SLIDE commandline incremental small

    $ ls 110217_SN167_0251_A8158LABXX/1_110217_A8158LABXX_barcode

    1_110217_A8158LABXX_5_1_fastq.txt
    1_110217_A8158LABXX_5_2_fastq.txt
    1_110217_A8158LABXX_6_1_fastq.txt
    1_110217_A8158LABXX_6_2_fastq.txt
    1_110217_A8158LABXX_unmatched_1_fastq.txt
    1_110217_A8158LABXX_unmatched_2_fastq.txt
    1_110217_A8158LABXX_bc.metrics
    SampleSheet-barcodes.cfg

    $ cat 1_110217_A8158LABXX_barcode/SampleSheet-barcodes.cfg 
    6 TGACCAA
    5 GCCAATA

    $ cat 1_110217_A8158LABXX_barcode/1_110217_A8158LABXX_bc.metrics
    5   25873823
    6   29331527
    unmatched   3051668
