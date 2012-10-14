import os
import sys
import subprocess
import unittest
import collections

from Bio import SeqIO

from scilifelab.utils.misc import safe_makedir

class SciLifeTest(unittest.TestCase):
    def setUp(self):
        self._install_1000g_test_files(os.path.join(os.path.dirname(__file__), "data", "production"))
        #self._install_bcbio_test_files(os.path.join(os.path.dirname(__file__), "data", "production"))

    def _install_1000g_test_files(self, data_dir):
        """Download 1000 genomes exome data

        See ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/sequence_indices/20120522.sequence.index
        for an index of recent sequencing runs. Load data into R

        df <- read.table("20120522.sequence.index", header=TRUE, sep="\t", fill=TRUE, as.is=TRUE)
        df$INDIVIDUAL =gsub("/sequence.*", "", gsub("data/", "", df$FASTQ_FILE))
        
        and select individual based on sequencing platform

        tapply(df$INSTRUMENT_MODEL, df$INDIVIDUAL, function(x) {levels(as.factor(x))})

        Here sequencing data from individual NA21137 has (arbitrarily) been chosen for download. Sequencing was
        done at BROAD institute on a Illumina HiSeq 2000.

        Method: download the fastq file and divide sequences into batches of 10000 
        sequences to emulate different projects. Then run data through pipeline to 
        generate downstream data. Downloading partial bam files with curl is possible, but
        then SamToFastq complains about unpaired mates so in the end it is better
        to download the entire bam (500M).

        The following code extracts 200000 reads that map to the 1.2 first Mb of 
        chromosome 11.
        """
        individual = "NA21137"
        base_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/{}".format(individual)
        bam_url = os.path.join(base_url, "exome_alignment", "NA21137.chrom11.ILLUMINA.bwa.GIH.exome.20111114.bam")
        tmp_dir = os.path.join(os.path.dirname(__file__), "tmp")
        if not os.path.exists(tmp_dir):
            safe_makedir(tmp_dir)
        bamfile = os.path.join(tmp_dir, os.path.basename(bam_url))
        ## Here check should be done on input files to pipeline; if not present, then
        ## download bamfile and generate fastq files
        if not os.path.exists(bamfile):
            cl = ["curl", bam_url, "-o", bamfile]
            subprocess.check_call(cl)

        ## Generate fastq files from bam
        if not os.getenv("PICARD_HOME", None):
            print "No environment variable PICARD_HOME set; exiting"
            sys.exit()
        try:
            cl = ["java", "-Xmx2g", "-XX:-UseGCOverheadLimit", "-jar", os.path.join(os.getenv("PICARD_HOME", os.curdir), "SamToFastq.jar"),
                  "INPUT={}".format(bamfile), "INCLUDE_NON_PF_READS=False", "FASTQ=reads_1.fq", "SECOND_END_FASTQ=reads_2.fq", "VALIDATION_STRINGENCY=SILENT"]
            if not os.path.exists("reads_1.fq"):
                subprocess.check_call(cl)
        except:
            print "Failed: {}".format(cl)
            raise

        ## Unfortunately the fastq files are not "paired". Here we
        ## loop the files and write to outfiles only if there are
        ## paired reads
        seqs1 = self._read_fastq("reads_1.fq")
        seqs2 = self._read_fastq("reads_2.fq")
        print "Writing fastq file 1..."
        self._write_fastq("seqs_1.fastq", seqs1, [x.id[0:-2] for x in seqs2])
        print "Writing fastq file 2..."
        self._write_fastq("seqs_2.fastq", seqs2, [x.id[0:-2] for x in seqs1])
        print "Done writing files"

    def _read_fastq(self, fn, numreads=250000):
        fh = open(fn, "rU")
        i = 0
        seqs = []
        for rec in SeqIO.parse(fh, "fastq"):
            i = i + 1
            if i % 10000==0:
                print "Read {} sequences...".format(i)
            ## For reads without a mate. SamToFastq should exclude
            ## these but apparently that doesn't happen
            if not rec.id[-2] == "/":
                print "excluding id {}".format(rec.id)
                continue
            seqs.append(rec)
            if i >= numreads:
                return seqs
    
    def _write_fastq(self, fn, seqs, ids):
        fh = open(fn, "w")
        i = 0
        for rec in seqs:
            if rec.id[-2] != "/":
                continue
            if rec.id[0:-2] in ids:
                del ids[ids.index(rec.id[0:-2])]
                i = i + 1
                SeqIO.write(rec, fh, "fastq")
                if i % 10000==0:
                    print "Wrote {} sequences...".format(i)
        fh.close()
            

    
