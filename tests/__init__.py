import os
import sys
import subprocess
import unittest
import logbook

from Bio import SeqIO

from scilifelab.utils.misc import safe_makedir

LOG = logbook.Logger(__name__)

## Directories and constants
filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
FLOWCELL = "120924_SN0002_0003_CC003CCCXX"
ARCHIVE = os.path.join(filedir, "data", "archive")
PRODUCTION = os.path.join(filedir, "data", "production")
GENOMES = os.path.join(filedir, "data", "genomes")

## Samplesheets - add when needed
SAMPLESHEETS = {}
SAMPLESHEETS['C003CCCXX'] = """FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
C003CCCXX,1,P001_101_index3,hg19,TGACCA,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,1,P001_102_index6,hg19,ACAGTG,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,2,P002_101_index3,hg19,TGACCA,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P002_102_index6,hg19,ACAGTG,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P002_103_index8,hg19,TGGTCA,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P003_101_index1,hg19,AGTGCG,J__Doe_00_03,N,R1,NN,J__Doe_00_03
C003CCCXX,2,P003_102_index2,hg19,TGTGCG,J__Doe_00_03,N,R1,NN,J__Doe_00_03
C003CCCXX,2,P003_103_index6,hg19,CGTTAA,J__Doe_00_03,N,R1,NN,J__Doe_00_03"""


def setUpModule():
    LOG.info("Running the setUpModule")
    _check_requirements()
    ## Add function to check existence of output files
    _install_1000g_test_files(os.path.join(os.path.dirname(__file__), "data", "production"))

def _check_requirements():
    """Assert the existence of paths and programs"""
    envvars = ["PICARD_HOME", "GATK_HOME"]
    for ev in envvars:
        if not os.getenv(ev, None):
            LOG.warn("Required environment variable {} missing".format({}))
            sys.exit()

def _install_1000g_test_files(data_dir):
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
    bam_url = os.path.join(base_url, "exome_alignment", "{}.chrom11.ILLUMINA.bwa.GIH.exome.20111114.bam".format(individual))
    tmpdir = os.path.join(os.path.dirname(__file__), "tmp")
    if not os.path.exists(tmpdir):
        safe_makedir(tmpdir)
    bamfile = os.path.join(tmpdir, os.path.basename(bam_url))
    ## Here check should be done on input files to pipeline; if not present, then
    ## download bamfile and generate fastq files
    if not os.path.exists(bamfile):
        cl = ["curl", bam_url, "-o", bamfile]
        subprocess.check_call(cl)

    _bam_to_fastq(bamfile, os.path.join(tmpdir, "reads"))
    r1 = os.path.join(tmpdir, "reads_1.fq")
    r2 = os.path.join(tmpdir, "reads_2.fq")
    _pair_fastq_files(r1, r2, os.path.join(tmpdir, "seqs"))
    _make_casava_archive_files()

def _bam_to_fastq(bamfile, out_prefix):
    """Convert bam to fastq file. Outputs paired reads"""
    try:
        cl = ["java", "-Xmx2g", "-XX:-UseGCOverheadLimit", "-jar", os.path.join(os.getenv("PICARD_HOME", os.curdir), "SamToFastq.jar"),
              "INPUT={}".format(bamfile), "INCLUDE_NON_PF_READS=False", "FASTQ={}_1.fq".format(out_prefix), "SECOND_END_FASTQ={}_2.fq".format(out_prefix), "VALIDATION_STRINGENCY=SILENT"]
        if not os.path.exists("{}_1.fq".format(out_prefix)):
            subprocess.check_call(cl)
    except:
        LOG.warn("Failed to run SamToFastq: {}".format(cl))
        raise

def _pair_fastq_files(r1, r2, out_prefix):
    """Pair fastq files. Unfortunately the fastq files are not
    "paired". Here we loop the files and write to outfiles only if
    there are paired reads
    """
    seqs1 = _read_fastq(r1)
    seqs2 = _read_fastq(r2)
    LOG.info("Writing fastq file 1...")
    _write_fastq("{}_1.fastq".format(out_prefix), seqs1, [x.id[0:-2] for x in seqs2])
    LOG.info("Writing fastq file 2...")
    _write_fastq("{}_2.fastq".format(out_prefix), seqs2, [x.id[0:-2] for x in seqs1])
    LOG.info("Done writing files")

def _read_fastq(fn, numreads=250000):
    fh = open(fn, "rU")
    i = 0
    seqs = []
    for rec in SeqIO.parse(fh, "fastq"):
        i = i + 1
        if i % 10000==0:
            LOG.info("Read {} sequences...".format(i))
        ## For reads without a mate. SamToFastq should exclude
        ## these but apparently that doesn't happen
        if not rec.id[-2] == "/":
            LOG.warning("excluding id {}".format(rec.id))
            continue
        seqs.append(rec)
        if i >= numreads:
            return seqs

def _write_fastq(fn, seqs, ids):
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
                LOG.info("Wrote {} sequences...".format(i))
    fh.close()

def _make_casava_archive_files(fc, ssname):
    fc_dir = os.path.join(ARCHIVE, fc)
    if not os.path.exists(fc_dir):
        safe_makedir(fc_dir)
    with open(os.path.join(fc_dir, "{}.csv".format(ssname)), "w") as fh:
        fh.write(SAMPLESHEETS[ssname])

def _download_genome_and_index(species="hg19", chr="11", start=None, end=None):
    """Download chromosome from ucsc, extract a given region and
    format for bwa and bowtie2.
    """
    pass

def _index_bwa(fn, outdir=None):
    """Index bwa"""
    pass
