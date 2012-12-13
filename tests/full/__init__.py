import os
import sys
import re
import subprocess
import unittest
import logbook
import gzip
from mako.template import Template
from cStringIO import StringIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

from scilifelab.utils.misc import safe_makedir, filtered_walk
from scilifelab.bcbio.flowcell import Flowcell
from scilifelab.pm.bcbio.utils import fc_id, fc_parts, fc_fullname

LOG = logbook.Logger(__name__)

## Directories and constants
filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
FLOWCELL = {"C003CCCXX": "120924_SN0002_0003_AC003CCCXX",
            "B002BBBXX": "121015_SN0001_0002_BB002BBBXX"}
ARCHIVE = os.path.join(filedir, "data", "archive")
PRODUCTION = os.path.join(filedir, "data", "production")
GENOMES = os.path.join(filedir, "data", "genomes")
CONFIG = os.path.join(filedir, "data", "config")
tmpdir = os.path.join(os.path.dirname(__file__), "tmp")
CURLFILESIZE = 30000000
NUMREADS = 50000

## Postprocess file
PPTEMPLATE = Template(filename=os.path.join(CONFIG, "post_process.mako"))
POSTPROCESS = os.path.join(CONFIG, "post_process.yaml")

## Samplesheets - add when needed
SAMPLESHEETS = {}
SAMPLESHEETS['C003CCCXX'] = """FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
C003CCCXX,1,P001_101_index3,hg19,TGACCA,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,2,P001_102_index6,hg19,ACAGTG,J__Doe_00_01,N,R1,NN,J__Doe_00_01
C003CCCXX,3,P001_101_index3,hg19,TGACCA,J__Doe_00_02,N,R1,NN_failed,J__Doe_00_02
C003CCCXX,3,P002_102_index6,hg19,ACAGTG,J__Doe_00_02,N,R1,NN,J__Doe_00_02
C003CCCXX,4,P003_101_index6,hg19,CGTTAA,J__Doe_00_03,N,R1,NN_failed,J__Doe_00_03"""
SAMPLESHEETS['B002BBBXX'] = """FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
B002BBBXX,1,P001_101_index3,hg19,TGACCA,J__Doe_00_01,N,R1,NN,J__Doe_00_01
B003BBBXX,2,P002_101_index3,hg19,TGACCA,J__Doe_00_02,N,R1,NN,J__Doe_00_02"""
NSAMPLES = sum([(len(SAMPLESHEETS[k].split("\n"))-1) for k in SAMPLESHEETS.keys()])
PROJECTS = ["J.Doe_00_01", "J.Doe_00_02", "J.Doe_00_03"]
RUNINFO=Template("""<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
  <Run Id="${flowcell}" Number="1">
    <Flowcell>${fc_id}</Flowcell>
    <Instrument>${instrument}</Instrument>
    <Date>${date}</Date>
    <Reads>
      <Read Number="1" NumCycles="101" IsIndexedRead="N" />
      <Read Number="2" NumCycles="7" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="101" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="8" SurfaceCount="2" SwathCount="3" TileCount="16" />
    <AlignToPhiX>
      <Lane>1</Lane>
      <Lane>2</Lane>
      <Lane>3</Lane>
      <Lane>4</Lane>
      <Lane>5</Lane>
      <Lane>6</Lane>
      <Lane>7</Lane>
      <Lane>8</Lane>
    </AlignToPhiX>
  </Run>
</RunInfo>
""")

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
    """Set up test files for scilifelab pipeline tests. The setup
    covers some typical situations, such as multiplexing, samples run
    on several flowcells, and same sample being run on several lanes
    in one flowcell.

    In short, the setup
    - downloads data from 1000 genomes (exome data from chr11, 0-2Mb)
    - generates fastq files in an archive folder
    - installs genome references (phix, hg19)
    - downloads dbsnp data for chr11, 0-2Mb
    - runs run_bcbb_pipeline.py -s to install fastq files to production folder
    - runs automated_initial_analysis.py
    """
    pattern = "14_write_metrics.txt"
    def filter_fn(f):
        return re.search(pattern, f) != None

    n = sum([len(filtered_walk(os.path.join(PRODUCTION, x), filter_fn)) for x in PROJECTS])
    if n == NSAMPLES:
        LOG.info("All samples have been run, requirements for downstream tests satisfied")
        return
    LOG.info("Running setUpModule")
    _check_requirements()
    ## Add function to check existence of output files
    _install_1000g_test_files(os.path.join(os.path.dirname(__file__), "data", "production"))

    _install_phix()
    dbsnp = _install_dbsnp_entrez()
    (omni_out, hapmap_out, mills_out) = _install_training_data()

    _download_ucsc_genome_and_index()
    ## Install post_process file
    fh = open(POSTPROCESS, "w")
    fh.write(PPTEMPLATE.render(**{'store_dir':ARCHIVE, 'base_dir':PRODUCTION, 'dbsnp':dbsnp, 'omni':omni_out, 'hapmap':hapmap_out, 'mills':mills_out}))
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

    ## Install files in production with run_bcbb_pipeline.py
    for k in FLOWCELL.keys():
        install = False
        for ss in SAMPLESHEETS[k].split("\n"):
            vals = ss.split(",")
            if vals[0]=="FCID":
                continue
            outdir = os.path.join(PRODUCTION, "{}".format(vals[5].replace("__", ".")), "{}".format(vals[2]), "{}_{}".format(FLOWCELL[k].split("_")[0],FLOWCELL[k].split("_")[-1]))
            r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
            r2 = os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1]))
            LOG.info("Looking for {} and {}".format(r1, r2))
            if not os.path.exists(r1) or not os.path.exists(r2):
                install = True
                break
        if install:
            LOG.info("Installing files with run_bcbb_pipeline.py for flowcell {}".format(k))
            cl = ["run_bcbb_pipeline.py", "-s", "-g", POSTPROCESS, os.path.join(ARCHIVE, FLOWCELL[k])]
            subprocess.check_call(cl)
        else:
            LOG.info("All files present; not running run_bcbb_pipeline.py")
    
    ## Run pipeline on samples 
    pattern = "-bcbb-config.yaml$"
    yamlfiles = []
    ## http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    ## [item for sublist in l for item in sublist]
    yamlfiles = [item for sublist in [filtered_walk(os.path.join(PRODUCTION, x), filter_fn) for x in PROJECTS] for item in sublist]
    orig_dir = os.path.abspath(os.curdir)
    for yamlconfig in yamlfiles:
        try:
            LOG.info("cding to {}".format(os.path.abspath(os.curdir)))
            os.chdir(os.path.dirname(yamlconfig))
            cl = ["automated_initial_analysis.py", POSTPROCESS, os.path.join(os.path.pardir, os.path.basename(os.path.dirname(yamlconfig))), yamlconfig]
            if not os.path.exists(os.path.join(os.path.dirname(yamlconfig), "14_write_metrics.txt")):
                LOG.info("Running pipeline: {}".format(" ".join(cl)))
                subprocess.check_call(cl)
        finally:
            os.chdir(orig_dir)
            LOG.info("Finished pipeline run and cd back to {}".format(orig_dir))


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
    ## FIXME: checks should be done on output from _pair_fastq_files
    smallbamfile = bamfile.replace(".bam", ".small.bam")
    if not os.path.exists(smallbamfile):
        LOG.info("downloading {} from {}".format(bamfile, base_url))
        cl = ["curl", bam_url, "-o", smallbamfile, "-r", "0-{}".format(CURLFILESIZE)]
        subprocess.check_call(cl)
        LOG.info("finished creating {}".format(smallbamfile))

    _bam_to_fastq(smallbamfile, os.path.join(tmpdir, "reads"))
    r1 = os.path.join(tmpdir, "reads_1.fq")
    r2 = os.path.join(tmpdir, "reads_2.fq")
    _pair_fastq_files(r1, r2, os.path.join(tmpdir, "seqs"))
    ## Install archive files
    _make_casava_archive_files(FLOWCELL["C003CCCXX"], "C003CCCXX", os.path.join(tmpdir, "seqs"))
    ## FIXME: startiter doesn't work, now generating identical files
    _make_casava_archive_files(FLOWCELL["B002BBBXX"], "B002BBBXX", os.path.join(tmpdir, "seqs"))

def _bam_to_fastq(bamfile, out_prefix):
    """Convert bam to fastq file. Outputs paired reads"""
    LOG.info("Converting bam {} to fastq with prefix {}".format(bamfile, out_prefix))
    try:
        cl = ["java", "-Xmx2g", "-XX:-UseGCOverheadLimit", "-jar", os.path.join(os.getenv("PICARD_HOME", os.curdir), "SamToFastq.jar"),
              "INPUT={}".format(bamfile), "INCLUDE_NON_PF_READS=False", "F={}_1.fq".format(out_prefix), "F2={}_2.fq".format(out_prefix), "VALIDATION_STRINGENCY=SILENT", "MAX_RECORDS_IN_RAM=5000000"]
        if not os.path.exists("{}_1.fq".format(out_prefix)):
            subprocess.check_call(cl)
    except:
        LOG.warn("Failed to run SamToFastq: {}".format(cl))
        LOG.info("This is expected since the input bamfile is truncated")
        pass

def _pair_fastq_files(r1, r2, out_prefix):
    """Pair fastq files. Unfortunately the fastq files are not
    "paired". Here we loop the files and write to outfiles only if
    there are paired reads
    """
    if os.path.exists("{}_1.fastq".format(out_prefix)):
        return
    LOG.info("reading {} and {}".format(r1, r2))
    seqs1 = _read_fastq(r1)
    seqs2 = _read_fastq(r2)
    LOG.info("Read {} sequences from read1, {} sequences from read2".format(len(seqs1), len(seqs1)))
    LOG.info("Writing fastq file 1...")
    _write_fastq("{}_1.fastq".format(out_prefix), seqs1, [x.id[0:-2] for x in seqs2])
    LOG.info("Writing fastq file 2...")
    _write_fastq("{}_2.fastq".format(out_prefix), seqs2, [x.id[0:-2] for x in seqs1])
    LOG.info("Done writing files")

def _read_fastq(fn, numreads=NUMREADS):
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
            i = i + 1
            SeqIO.write(rec, fh, "fastq")
            if i % 10000==0:
                LOG.info("Wrote {} sequences...".format(i))
    fh.close()

def _make_casava_archive_files(fc, ssname, prefix, startiter = 1, nseqout=1000):
    fc_dir = os.path.join(ARCHIVE, fc)
    if not os.path.exists(fc_dir):
        safe_makedir(fc_dir)
    with open(os.path.join(fc_dir, "{}.csv".format(ssname)), "w") as fh:
        fh.write(SAMPLESHEETS[ssname])
    with open(os.path.join(fc_dir, "RunInfo.xml"), "w") as fh:
        fh.write(RUNINFO.render(**{'flowcell':os.path.basename(fc), 'fc_id':fc_id(fc), 'date':fc_parts(fc)[0], 'instrument':split("_", fc)[1]}))
    outf1 = []
    outf2 = []
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
        if len(vals) == 0:
            continue
        outdir = os.path.join(fc_dir, "Unaligned", "Project_{}".format(vals[5]), "Sample_{}".format(vals[2]))
        if not os.path.exists(outdir):
            safe_makedir(outdir)
        with open(os.path.join(outdir, "SampleSheet.csv"), "w") as fh:
            LOG.info("Writing to {}".format(os.path.join(outdir, "SampleSheet.csv")))
            fh.write("{}\n".format(header))
            fh.write("{}\n".format(row))
        r1 = os.path.join(outdir, "{}_{}_L00{}_R1_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        r2 = os.path.join(outdir, "{}_{}_L00{}_R2_001.fastq.gz".format(vals[2], vals[4], vals[1]))
        if os.path.exists(r1):
            LOG.info("{} already exists: if you want to rerun file generation remove {}".format(r1, r1))
            return 
        outf1.append(r1)
        outf2.append(r2)

    ## Write sequences
    with open("{}_1.fastq".format(prefix), "r") as fh:
        _write_sample_fastq(fh, outf1, startiter=startiter, nseqout=nseqout)
    with open("{}_2.fastq".format(prefix), "r") as fh:
        _write_sample_fastq(fh, outf2, startiter=startiter, nseqout=nseqout)

def _write_sample_fastq(fh, outfiles, startiter=0, nseqout=1000):
    i = 0
    j = 0 - startiter
    outh = []
    for of in outfiles:
        LOG.info("Opening gzip file {}".format(of))
        oh = gzip.open(of, "w")
        outh.append(oh)
    n = len(outh)
    totseqout = n * nseqout
    LOG.info("writing {} sequences per file for {} samples, {} sequences in total".format(nseqout, n, totseqout))
    for rec in SeqIO.parse(fh, "fastq"):
        j = j + 1
        if j < 0:
            continue
        SeqIO.write(rec, outh[i%n], "fastq")
        i = i + 1
        if (i % 1000 == 0):
            LOG.info("read {} sequences from {}".format(i, fh.name))
        if i > totseqout:
            break
    [h.close() for h in outh]

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
        LOG.info("Downloading {} from {} with curl".format(os.path.join(genomedir, os.path.basename(url)), url))
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

def _index_bwa(fn, label="bwa"):
    """Index bwa"""
    LOG.info("Indexing {} with bwa".format(fn))
    outdir = os.path.join(os.path.dirname(fn), os.pardir,label)
    if not os.path.exists(outdir):
        safe_makedir(outdir)
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.amb".format( os.path.basename(fn)))):
        LOG.info("{} exists; not doing anything".format(fn))
        return os.path.join(outdir, os.path.basename(fn))
    cl = ["bwa", "index", os.path.abspath(os.path.join(outdir, os.path.basename(fn)))]
    subprocess.check_call(cl)
    LOG.info("Finished indexing {} with bwa".format(fn))
    return os.path.join(outdir, os.path.basename(fn))

def _index_bowtie(fn, label="bowtie"):
    """Index bowtie"""
    LOG.info("Indexing {} with bowtie".format(fn))
    outdir = os.path.join(os.path.dirname(fn), os.pardir, label)
    if not os.path.exists(outdir):
        safe_makedir(outdir)
        os.symlink(fn, os.path.join(outdir, os.path.basename(fn)))
    if os.path.exists(os.path.join(outdir,"{}.1.ebwt".format( os.path.splitext(os.path.basename(fn))[0]))):
        LOG.info("{} exists; not doing anything".format(fn))
        return os.path.splitext(os.path.join(outdir, os.path.basename(fn)))[0]
    cl = ["bowtie-build", os.path.abspath(os.path.join(outdir, os.path.basename(fn))), os.path.splitext(os.path.abspath(os.path.join(outdir, os.path.basename(fn))))[0]]
    subprocess.check_call(cl)
    LOG.info("Finished indexing {} with bowtie".format(fn))
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
    LOG.info("Installing phix")
    build = "phix"
    genomedir = os.path.join(GENOMES, genomes[build]['species'], build, "seq")
    fn = os.path.join(genomedir, "phix.fa")
    if not os.path.exists(genomedir):
        LOG.info("Creating {}".format(genomedir))
        safe_makedir(genomedir)
    if not os.path.exists(fn):
        try:
            LOG.info("Opening file {}".format(fn))
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

##############################
## Variation data
##############################
dbsnp_header = """##fileformat=VCFv4.0
##fileDate=20101103
##source=dbSNP
##dbSNP_BUILD_ID=132
##reference=GRCh37
##phasing=partial
##variationPropertyDocumentationUrl=ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf
##INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=VP,Number=1,Type=String,Description="Variation Property">
##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description="First SNP Build for RS">
##INFO=<ID=WGT,Number=1,Type=Integer,Description="Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more">
##INFO=<ID=VC,Number=1,Type=String,Description="Variation Class">
##INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
##INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
##INFO=<ID=TPA,Number=0,Type=Flag,Description="Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)">
##INFO=<ID=PMC,Number=0,Type=Flag,Description="Links exist to PubMed Central article">
##INFO=<ID=S3D,Number=0,Type=Flag,Description="Has 3D structure - SNP3D table">
##INFO=<ID=SLO,Number=0,Type=Flag,Description="Has SubmitterLinkOut - From SNP->SubSNP->Batch.link_out">
##INFO=<ID=NSF,Number=0,Type=Flag,Description="Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44">
##INFO=<ID=NSM,Number=0,Type=Flag,Description="Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42">
##INFO=<ID=NSN,Number=0,Type=Flag,Description="Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41">
##INFO=<ID=REF,Number=0,Type=Flag,Description="Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8">
##INFO=<ID=SYN,Number=0,Type=Flag,Description="Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3">
##INFO=<ID=U3,Number=0,Type=Flag,Description="In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53">
##INFO=<ID=U5,Number=0,Type=Flag,Description="In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55">
##INFO=<ID=ASS,Number=0,Type=Flag,Description="In acceptor splice site FxnCode = 73">
##INFO=<ID=DSS,Number=0,Type=Flag,Description="In donor splice-site FxnCode = 75">
##INFO=<ID=INT,Number=0,Type=Flag,Description="In Intron FxnCode = 6">
##INFO=<ID=R3,Number=0,Type=Flag,Description="In 3' gene region FxnCode = 13">
##INFO=<ID=R5,Number=0,Type=Flag,Description="In 5' gene region FxnCode = 15">
##INFO=<ID=OTH,Number=0,Type=Flag,Description="Has other snp with exactly the same set of mapped positions on NCBI refernce assembly.">
##INFO=<ID=CFL,Number=0,Type=Flag,Description="Has Assembly conflict. This is for weight 1 and 2 snp that maps to different chromosomes on different assemblies.">
##INFO=<ID=ASP,Number=0,Type=Flag,Description="Is Assembly specific. This is set if the snp only maps to one assembly">
##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
##INFO=<ID=VLD,Number=0,Type=Flag,Description="Is Validated.  This bit is set if the snp has 2+ minor allele count based on frequency or genotype data.">
##INFO=<ID=G5A,Number=0,Type=Flag,Description=">5% minor allele frequency in each and all populations">
##INFO=<ID=G5,Number=0,Type=Flag,Description=">5% minor allele frequency in 1+ populations">
##INFO=<ID=HD,Number=0,Type=Flag,Description="Marker is on high density genotyping kit (50K density or greater).  The snp may have phenotype associations present in dbGaP.">
##INFO=<ID=GNO,Number=0,Type=Flag,Description="Genotypes available. The snp has individual genotype (in SubInd table).">
##INFO=<ID=KGPilot1,Number=0,Type=Flag,Description="1000 Genome discovery(pilot1) 2009">
##INFO=<ID=KGPilot123,Number=0,Type=Flag,Description="1000 Genome discovery all pilots 2010(1,2,3)">
##INFO=<ID=KGVAL,Number=0,Type=Flag,Description="1000 Genome validated by second method">
##INFO=<ID=KGPROD,Number=0,Type=Flag,Description="1000 Genome production phase">
##INFO=<ID=PH1,Number=0,Type=Flag,Description="Phase 1 genotyped: filtered, non-redundant">
##INFO=<ID=PH2,Number=0,Type=Flag,Description="Phase 2 genotyped: filtered, non-redundant">
##INFO=<ID=PH3,Number=0,Type=Flag,Description="Phase 3 genotyped: filtered, non-redundant">
##INFO=<ID=CDA,Number=0,Type=Flag,Description="Variation is interrogated in a clinical diagnostic assay">
##INFO=<ID=LSD,Number=0,Type=Flag,Description="Submitted from a locus-specific database">
##INFO=<ID=MTP,Number=0,Type=Flag,Description="Microattribution/third-party annotation(TPA:GWAS,PAGE)">
##INFO=<ID=OM,Number=0,Type=Flag,Description="Has OMIM/OMIA">
##INFO=<ID=NOC,Number=0,Type=Flag,Description="Contig allele not present in SNP allele list. The reference sequence allele at the mapped position is not present in the SNP allele list, adjusted for orientation.">
##INFO=<ID=WTD,Number=0,Type=Flag,Description="Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory">
##INFO=<ID=NOV,Number=0,Type=Flag,Description="Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common.">
##INFO=<ID=GCF,Number=0,Type=Flag,Description="Has Genotype Conflict Same (rs, ind), different genotype.  N/N is not included.">
"""

omni="""##fileformat=VCFv4.1
##FILTER=<ID=NOT_POLY_IN_1000G,Description="Alternate allele count = 0">
##FILTER=<ID=badAssayMapping,Description="The mapping information for the SNP assay is internally inconsistent in the chip metadata">
##FILTER=<ID=dup,Description="Duplicate assay at same position with worse Gentrain Score">
##FILTER=<ID=id10,Description="Within 10 bp of an known indel">
##FILTER=<ID=id20,Description="Within 20 bp of an known indel">
##FILTER=<ID=id5,Description="Within 5 bp of an known indel">
##FILTER=<ID=id50,Description="Within 50 bp of an known indel">
##FILTER=<ID=refN,Description="Reference base is N. Assay is designed for 2 alt alleles">
##FORMAT=<ID=GC,Number=.,Type=Float,Description="Gencall Score">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[/gap/birdsuite/1kg/0.928975161471502.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false enable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=CR,Number=.,Type=Float,Description="SNP Callrate">
##INFO=<ID=GentrainScore,Number=.,Type=Float,Description="Gentrain Score">
##INFO=<ID=HW,Number=.,Type=Float,Description="Hardy-Weinberg Equilibrium">
##reference=human_g1k_v37.fasta
##source=infiniumFinalReportConverterV1.0
"""

hapmap = """##fileformat=VCFv4.1
##CombineVariants="analysis_type=CombineVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta rodBind=[/broad/shptmp/0.516962905488075.ASW.vcf, /broad/shptmp/0.516962905488075.CEU.vcf, /broad/shptmp/0.516962905488075.CHB.vcf, /broad/shptmp/0.516962905488075.CHD.vcf, /broad/shptmp/0.516962905488075.GIH.vcf, /broad/shptmp/0.516962905488075.JPT.vcf, /broad/shptmp/0.516962905488075.LWK.vcf, /broad/shptmp/0.516962905488075.MEX.vcf, /broad/shptmp/0.516962905488075.MKK.vcf, /broad/shptmp/0.516962905488075.TSI.vcf, /broad/shptmp/0.516962905488075.YRI.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub genotypemergeoption=UNSORTED variantmergeoption=UNION rod_priority_list=ASW,YRI,LWK,CHD,CHB,CEU,GIH,MKK,MEX,JPT,TSI printComplexMerges=false filteredAreUncalled=false minimalVCF=false setKey=set"
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[/broad/shptmp/ebanks//0.764768180511059.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##VariantsToVCF="analysis_type=VariantsToVCF input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=[X] excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta rodBind=[/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/rawdata/genotypes_chrX_ASW_phase3.3_consensus.b36_fwd.txt] rodToIntervalTrackName=null BTI_merge_rule=UNION DBSNP=/humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=1.0E-4 useOriginalQualities=false validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null logging_level=INFO log_to_file=null quiet_output_mode=false debug_mode=false help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sample=null"
##reference=Homo_sapiens_assembly18.fasta
##source=VariantsToVCF
"""
mills = """##fileformat=VCFv4.1
##CombineVariants="analysis_type=CombineVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[./indel_hg19_051711_leftAligned.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub genotypemergeoption=PRIORITIZE filteredrecordsmergetype=KEEP_IF_ANY_UNFILTERED rod_priority_list=mills printComplexMerges=false filteredAreUncalled=false minimalVCF=false setKey=set assumeIdenticalSamples=false minimumN=1 masterMerge=false mergeInfoWithMaxAC=false"
##FilterLiftedVariants="analysis_type=FilterLiftedVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta rodBind=[/broad/shptmp/delangel/tmp/0.101786952306615.sorted.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##LeftAlignVariants="analysis_type=LeftAlignVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_b36_both.fasta rodBind=[./indel_hg18_051711.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub"
##ValidateVariants="analysis_type=ValidateVariants input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_g1k_v37.fasta rodBind=[./indel_hg19_051711_leftAligned_collapsed.vcf] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub validationType=ALL doNotValidateFilteredRecords=false warnOnErrors=true"
##VariantsToVCF="analysis_type=VariantsToVCF input_file=[] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=null excludeIntervals=null reference_sequence=/humgen/1kg/reference/human_b36_both.fasta rodBind=[./indel_hg18_051711_sorted.txt] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=1 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sample=null fixRef=true"
"""

vcfheader="{}\n".format("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]))

def _install_dbsnp_entrez(build="hg19"):
    """Install a subset of snps using Entrez queries"""
    LOG.info("Installing dbsnp file for {}".format(genomes[build]['species']))
    variationdir = os.path.join(GENOMES, genomes[build]['species'], build, "variation")
    if not os.path.exists(variationdir):
        safe_makedir(variationdir)
    fn = os.path.join(variationdir, "dbsnp132_chr11.vcf")
    if not os.path.exists(fn):
        try:
            # http://www.ncbi.nlm.nih.gov/books/NBK44454/#Search.how_do_i_search_dbsnp_for_the_tot
            ## This will actually only download a subset of snps
            handle = Entrez.esearch(db="snp", retmax=8000, term="\"Homo sapiens\"[Organism] AND (11[CHR] AND (1[CHRPOS] : 2000000[CHRPOS])")
            record = Entrez.read(handle)
            records = []
            ## For some reason the first entries are more or less empty
            start = 4000
            delta = 200
            for i in xrange(start, len(record['IdList']), delta):
                LOG.info("retrieving dbsnp records {} - {}".format(i, i+delta))
                h = Entrez.efetch(db="snp", id=record['IdList'][i:i+delta], rettype="flt", retmax=delta)
                lbuffer = None
                lines = []
                while True:
                    if lbuffer:
                        l = lbuffer
                        lbuffer = None
                    else:
                        l = h.readline()
                    if l.startswith("rs") or len(lines) > 20:
                        if lines:
                            lbuffer = l
                            rec = _dbsnp_line(lines)
                            if rec is None:
                                break
                            records.append(rec)
                            lines = []
                        else:
                            lines.append(l.rstrip())
                    else:
                        lines.append(l.rstrip())
        except:
            LOG.warning("Entrez query failed")
            pass
        LOG.info("Writing file {}".format(fn))
        fh = open(fn, "w")
        fh.write(dbsnp_header)
        ## Header must be tab-separated, otherwise GATK complains...
        fh.write(vcfheader)
        for rec in sorted(records, key=lambda x: int(x.split()[1])):
            fh.write(rec)
            fh.write("\n")
        fh.close()

    return fn
    
def _dbsnp_line(lines):
    if not lines:
        return None
    line = "".join(lines)
    try:
        m = re.search("^(rs[0-9]+).*alleles=([A-Z]+)/([A-Z]+).*chr=([0-9]+).*chr-pos=([0-9]+)", line)
    except:
        LOG.warning("regexp search failed")
        raise
    if m:
        return "\t".join(["chr{}".format(m.groups()[3]), m.groups()[4], m.groups()[0], m.groups()[1], m.groups()[2], ".", ".", "dbSNPBuildID=132;VP=050000020005000000000100;WGT=1"])
    else:
        return None

def _install_dbsnp_file(build="hg19"):
    """Download a (large) dbsnp file and extract a region from chr 11"""
    variationdir = os.path.join(GENOMES, genomes[build]['species'], build, "variation")
    url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/dbsnp132_20101103.vcf.gz"
    fn = os.path.join(variationdir, "dbsnp132.vcf.gz")
    dbsnp = os.path.join(variationdir, "dbsnp132_chr11.vcf")
    if not os.path.exists(variationdir):
        safe_makedir(variationdir)
    try:
        LOG.info("Downloading {} from {} with curl".format(fn, url))
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

def _install_training_data(build="hg19"):
    variationdir = os.path.join(GENOMES, genomes[build]['species'], build, "variation")
    omni_out = os.path.join(variationdir, "1000G_omni2.5.vcf")
    hapmap_out = os.path.join(variationdir, "hapmap_3.3.vcf")
    mills_out = os.path.join(variationdir, "Mills_Devine_2hit.indels.vcf")
    fh = open(omni_out, "w")
    fh.write(omni)
    fh.write(vcfheader)
    fh.close()
    fh = open(hapmap_out, "w")
    fh.write(hapmap)
    fh.write(vcfheader)
    fh.close()
    fh = open(mills_out, "w")
    fh.write(mills)
    fh.write(vcfheader)
    fh.close()
    return (omni_out, hapmap_out, mills_out)
