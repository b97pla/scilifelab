import os
import subprocess
import unittest
import collections

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
            #cl = ["curl", bam_url, "-r", "0-20000000", "-o", bamfile]
            cl = ["curl", bam_url, "-o", bamfile]
            subprocess.check_call(cl)
        
        samfile = os.path.join(tmp_dir, os.path.basename(bam_url).replace(".bam", ".part.sam"))
        cl = ["samtools", "view", "-h", bamfile, "|", "head", "-200000", ">", samfile]
        subprocess.Popen(" ".join(cl), shell=True)
        bamfile2 = os.path.join(tmp_dir, os.path.basename(bam_url).replace(".bam", ".part.bam"))
        cl = ["samtools", "view", "-Sbh", samfile, ">", bamfile2]
        subprocess.Popen(" ".join(cl), shell=True)
        cl = ["samtools", "index", bamfile2]
        subprocess.check_call(cl)

    def _install_bcbio_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        
        DlInfo = collections.namedtuple("DlInfo", "fname dirname version")
        download_data = [DlInfo("110106_FC70BUKAAXX.tar.gz", None, None),
                         DlInfo("genomes_automated_test.tar.gz", "genomes", 6),
                         DlInfo("110907_ERP000591.tar.gz", None, None),
                         DlInfo("100326_FC6107FAAXX.tar.gz", None, 2)]
        for dl in download_data:
            url = "http://chapmanb.s3.amazonaws.com/{fname}".format(fname=dl.fname)
            dirname = os.path.join(data_dir, 
                                   dl.fname.replace(".tar.gz", "") if dl.dirname is None
                                   else dl.dirname)
            if os.path.exists(dirname) and dl.version is not None:
                version_file = os.path.join(dirname, "VERSION")
                is_old = True
                if os.path.exists(version_file):
                    with open(version_file) as in_handle:
                        version = int(in_handle.read())
                    is_old = version < dl.version
                if is_old:
                    shutil.rmtree(dirname)
            if not os.path.exists(dirname):
                self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print "Downloading data to {}".format(dirname)
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        if not os.path.exists(os.path.dirname(dirname)):
            safe_makedir(os.path.dirname(dirname))
        os.rename(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))
