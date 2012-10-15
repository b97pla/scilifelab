import os
import sys
import subprocess
import unittest
import logbook
import gzip
from mako.template import Template
from cStringIO import StringIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

from scilifelab.utils.misc import safe_makedir
from scilifelab.bcbio.flowcell import Flowcell

LOG = logbook.Logger(__name__)

## Directories and constants
filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
FLOWCELL = "120924_SN0002_0003_CC003CCCXX"
ARCHIVE = os.path.join(filedir, "data", "archive")
PRODUCTION = os.path.join(filedir, "data", "production")
GENOMES = os.path.join(filedir, "data", "genomes")
CONFIG = os.path.join(filedir, "data", "config")
tmpdir = os.path.join(os.path.dirname(__file__), "tmp")

## Postprocess file
PPTEMPLATE = Template(filename=os.path.join(CONFIG, "post_process.mako"))
POSTPROCESS = os.path.join(CONFIG, "post_process.yaml")

## Samplesheets - add when needed
SAMPLESHEETS = {}
SAMPLESHEETS['C003CCCXX'] = """FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
C003CCCXX,1,P001_101_index3,hg19,TGACCA,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,1,P001_102_index6,hg19,ACAGTG,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,2,P002_101_index3,hg19,TGACCA,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P002_102_index6,hg19,ACAGTG,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P002_103_index8,hg19,TGGTCA,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,2,P003_101_index1,hg19,AGTGCG,J__Doe_00_03,N,R1,NN,J__Doe_00_03
C003CCCXX,2,P003_101_index2,hg19,TGTGCG,J__Doe_00_03,N,R1,NN,J__Doe_00_03
C003CCCXX,2,P003_101_index6,hg19,CGTTAA,J__Doe_00_03,N,R1,NN,J__Doe_00_03"""

## Genome metadata
genomes = {'hg19':{'species':'Hsapiens', 'label':'Human (hg19)'},
           'phix':{'species':'phix', 'label':'phiX174'}}

## Index files
index_files = {'sam':{'file':os.path.join(CONFIG, "tool-data", "sam_fa_indices.loc"), 'data':StringIO()},
               'bwa':{'file':os.path.join(CONFIG, "tool-data", "bwa_index.loc"), 'data':StringIO()},
               'bowtie':{'file':os.path.join(CONFIG, "tool-data", "bowtie_indices.loc"), 'data':StringIO()},
               'bowtie2':{'file':os.path.join(CONFIG, "tool-data", "bowtie2_indices.loc"), 'data':StringIO()},
               'liftOver':{'file':os.path.join(CONFIG, "tool-data", "liftOver.loc"), 'data':StringIO()}
               }


def setUpModule():
    LOG.info("Running setUpModule")
    _check_requirements()
    ## Add function to check existence of output files
    _install_1000g_test_files(os.path.join(os.path.dirname(__file__), "data", "production"))
    _install_phix()
    dbsnp = _install_dbsnp()
    _download_ucsc_genome_and_index()
    ## Install post_process file
    fh = open(POSTPROCESS, "w")
    fh.write(PPTEMPLATE.render(**{'store_dir':ARCHIVE, 'base_dir':PRODUCTION, 'dbsnp':dbsnp}))
    fh.close()
    ## Install index files
    for k, v in index_files.iteritems():
        if not os.path.exists(os.path.dirname(v['file'])):
            safe_makedir(os.path.dirname(v['file']))
        fh = open(v['file'], "w")
        fh.write(v['data'].getvalue())
        fh.close()
    ## Make production dir
    if not os.path.exists(PRODUCTION):
        safe_makedir(PRODUCTION)


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
    _make_casava_archive_files(FLOWCELL, "C003CCCXX", os.path.join(tmpdir, "seqs"))

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
    if os.path.exists("{}_1.fastq".format(out_prefix)):
        return
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

def _make_casava_archive_files(fc, ssname, prefix, nseqout=100000):
    fc_dir = os.path.join(ARCHIVE, fc)
    if not os.path.exists(fc_dir):
        safe_makedir(fc_dir)
    with open(os.path.join(fc_dir, "{}.csv".format(ssname)), "w") as fh:
        fh.write(SAMPLESHEETS[ssname])
    #fc = Flowcell(os.path.join(fc_dir, "{}.csv".format(ssname)))
    h1 = open("{}_1.fastq".format(prefix), "r")
    h2 = open("{}_2.fastq".format(prefix), "r")
    outh1 = []
    outh2 = []
    basecall_stats_dir = os.path.join(fc_dir, "Unaligned", "Basecall_Stats_{}".format(ssname))
    if not os.path.exists(basecall_stats_dir):
        safe_makedir(basecall_stats_dir)
    for d in [os.path.join(basecall_stats_dir, x) for x in ["css", "Plots"]]:
        if not os.path.exists(d):
            safe_makedir(d)
    
    for row in SAMPLESHEETS[ssname].split("\n"):
        vals = row.split(",")
        if vals[0] == "FCID":
            header = row
            continue
        outdir = os.path.join(fc_dir, "Unaligned", "Project_{}".format(vals[5]), "Sample_{}".format(vals[2]))
        if not os.path.exists(outdir):
            safe_makedir(outdir)
        with open(os.path.join(outdir, "SampleSheet.csv"), "w") as fh:
            LOG.info("Writing to {}".format(os.path.join(outdir, "SampleSheet.csv")))
            fh.write("{}\n".format(header))
            fh.write("{}\n".format(row))
        r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        if os.path.exists(r1):
            LOG.info("{} already exists: if you want to rerun file generation remove {}".format(r1, r1))
            return 
        oh1 = gzip.open(r1, "w")
        outh1.append(oh1)
        oh2 = gzip.open(os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1])),"w")
        outh2.append(oh2)


    ## Write sequences
    i = 0
    n = len(outh1)
    for rec in SeqIO.parse(h1, "fastq"):
        SeqIO.write(rec, outh1[i % n], "fastq")
        i = i + 1
        if (i % 10000 == 0):
            LOG.info("read {} sequences from {}".format(i, h1.name))
        if i > nseqout:
            break
    [h.close() for h in outh1]
    i = 0
    n = len(outh2)
    for rec in SeqIO.parse(h2, "fastq"):
        SeqIO.write(rec, outh2[i%n], "fastq")
        i = i + 1
        if (i % 10000 == 0):
            LOG.info("read {} sequences from {}".format(i, h2.name))
        if i > nseqout:
            break

    [h.close() for h in outh2]
        
    h1.close()
    h2.close()

def _download_ucsc_genome_and_index(build="hg19", chr="chr11", start=0, end=2000000):
    """Download chromosome from ucsc, extract a given region and
    format for bwa and bowtie2.
    """
    ucsc = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/chromosomes/".format(build)
    url = os.path.join(ucsc, "{}.fa.gz".format(chr))
    genomedir = os.path.join(GENOMES, genomes[build]['species'], build, "seq")
    if not os.path.exists(genomedir):
        safe_makedir(genomedir)
    try:
        cl = ["curl", url, "-o", os.path.join(genomedir, os.path.basename(url))]
        if not os.path.exists(os.path.join(genomedir, os.path.basename(url))):
            subprocess.check_call(cl)
    except:
        pass
    if not os.path.exists(os.path.join(genomedir, os.path.basename(url).replace(".gz", ""))):
        rec = SeqIO.read(gzip.open(os.path.join(genomedir, os.path.basename(url)), "r"), "fasta")
        outh = open(os.path.join(genomedir, os.path.basename(url).replace(".gz", "")), "w")
        SeqIO.write(SeqRecord(rec.seq[start:end], rec.id, '', ''), outh, "fasta")
        outh.close()
    outfile = _index_bwa(os.path.join(genomedir, os.path.basename(url).replace(".gz", "")))
    index_files['sam']['data'].write("index\t{}\t{}\n".format(build, os.path.join(genomedir, os.path.basename(url).replace(".gz", ""))))
    index_files['bwa']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    outfile = _index_bowtie(os.path.join(genomedir, os.path.basename(url).replace(".gz", "")))
    index_files['bowtie']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    #outfile = _index_bowtie2(os.path.join(genomedir, os.path.basename(url).replace(".gz", "")))
    #index_files['bowtie2']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))

def _index_bwa(fn, label="bwa"):
    """Index bwa"""
    outdir = os.path.join(os.path.dirname(fn), os.pardir,label)
    if not os.path.exists(outdir):
        safe_makedir(outdir)
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.amb".format( os.path.basename(fn)))):
        return os.path.join(outdir, os.path.basename(fn))
    cl = ["bwa", "index", os.path.abspath(os.path.join(outdir, os.path.basename(fn)))]
    subprocess.check_call(cl)
    return os.path.join(outdir, os.path.basename(fn))

def _index_bowtie(fn, label="bowtie"):
    """Index bowtie"""
    outdir = os.path.join(os.path.dirname(fn), os.pardir, label)
    if not os.path.exists(outdir):
        safe_makedir(outdir)
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.1.ebwt".format( os.path.splitext(os.path.basename(fn))[0]))):
        return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]
    cl = ["bowtie-build", os.path.abspath(os.path.join(outdir, os.path.basename(fn))), os.path.splitext(os.path.abspath(os.path.join(outdir, os.path.basename(fn))))[0]]
    subprocess.check_call(cl)
    return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]

def _index_bowtie2(fn, label="bowtie2"):
    """Index bowtie2"""
    outdir = os.path.join(os.path.dirname(fn), os.pardir, label)
    if not os.path.exists(outdir):
        safe_makedir(outdir)
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.1.bt2".format( os.path.splitext(os.path.basename(fn))[0]))):
        return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]
    cl = ["bowtie2-build", os.path.abspath(os.path.join(outdir, os.path.basename(fn))), os.path.splitext(os.path.abspath(os.path.join(outdir, os.path.basename(fn))))[0]]
    subprocess.check_call(cl)
    return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]

def _install_phix():
    build = "phix"
    genomedir = os.path.join(GENOMES, genomes[build]['species'], build, "seq")
    fn = os.path.join(genomedir, "phix.fa")
    if not os.path.exists(genomedir):
        safe_makedir(genomedir)
    if not os.path.exists(fn):
        try:
            fh = open(fn, "w")
            handle = Entrez.efetch(db="nucleotide", id="9626372", rettype="fasta", retmode="text")
            rec = "".join(handle.readlines())
            fh.write(rec)
            fh.close()
        except:
            pass
    outfile = _index_bwa(fn, label="bwa")
    index_files['sam']['data'].write("index\t{}\t{}\n".format(build, fn))
    index_files['bwa']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    outfile = _index_bowtie(fn, label="bowtie")
    index_files['bowtie']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))
    #outfile = _index_bowtie2(fn, label="bowtie2")
    #index_files['bowtie2']['data'].write("{}\t{}\t{}\t{}\n".format(build, build, genomes[build]['label'], outfile))

def _install_dbsnp(build="hg19"):
    """Download a (large) dbsnp file and extract a region from chr 11"""
    variationdir = os.path.join(GENOMES, genomes[build]['species'], build, "variation")
    url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/dbsnp132_20101103.vcf.gz"
    fn = os.path.join(variationdir, "dbsnp132.vcf.gz")
    dbsnp = os.path.join(variationdir, "dbsnp132_chr11.vcf")
    if not os.path.exists(variationdir):
        safe_makedir(variationdir)
    try:
        cl = ["curl", url, "-o", fn]
        if not os.path.exists(os.path.join(variationdir, os.path.basename(fn))):
            subprocess.check_call(cl)
    except:
        pass
    fh = gzip.open(fn, "read")
    if not os.path.exists(dbsnp):
        of = open(dbsnp, "w")
        for r in fh:
            vals = r.split()
            if r.startswith("#"):
                of.write(r)
            if vals[0] != "11":
                continue
            if int(vals[1]) < 2000000:
                vals[0] = "chr{}".format(vals[0])
                of.write("\t".join(vals))
                of.write("\n")
            else:
                break
        of.close()
    return dbsnp
## run_bcbb_pipeline.py flowcell_path post_process.yaml [-custom runinfo -g]
