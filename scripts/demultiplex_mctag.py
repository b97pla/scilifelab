"""
Demultiplex haloplex data including molecular tags.
"""
# TODO add pairwise alignment of indexes to correct for sequencing error (biopython's Bio.pairwise2)
# TODO add directory parsing
# TODO support single-read runs (i.e. not paired-end)

from __future__ import print_function

import argparse
import Bio
import collections
import itertools
import os
import re
import sys

from scilifelab.illumina.hiseq import HiSeqSampleSheet
from scilifelab.utils.fastq_utils import FastQParser, FastQWriter

class NoIndexMatchException(Exception):
    def __init__(self, message, *args):
        if not message:
            message = "no match to supplied indexes."
        super(NoIndexMatchException, self).__init__(self, message)
        #Exception.__init__(self, message)
        self.args = args

class RevComMatchException(NoIndexMatchException):
    def __init__(self, message, *args):
        if not message:
            message = "no match to supplied indexes but match to reverse complement of supplied indexes. Double-check read direction."
        Exception.__init__(self, message)
        self.args = args



def main(read_one, read_two, read_index, sample_sheet, halo_index_file, halo_index_length, molecular_tag_length=None):

    halo_index_dict, halo_index_revcom_dict = load_indexes(halo_index_file)
    parse_readset(read_one, read_two, read_index, halo_index_dict, halo_index_revcom_dict, halo_index_length, molecular_tag_length)


def parse_readset(read_1_fq, read_2_fq, read_index_fq, halo_index_dict, halo_index_revcom_dict, halo_index_length, molecular_tag_length=None):
    """
    Parse input fastq files by index reads.
    """
    # TODO make a real directory, not this hack for testing
    output_dir = os.path.join(os.getcwd(), "processed_files")
    os.makedirs(output_dir)

    fqp_1, fqp_2, fqp_ind = map(FastQParser, (read_1_fq, read_2_fq, read_index_fq))

    match, non_match, revcom_match = 0, 0, 0

    for read_1, read_2, read_ind in itertools.izip(fqp_1, fqp_2, fqp_ind):
        # Get the haloplex index and the molecular tag from the index read's sequence data
        try:
            halo_index, molecular_tag = parse_index(read_ind[1], halo_index_dict, halo_index_revcom_dict,
                                                    halo_index_length, molecular_tag_length)
            match += 1
        except RevComMatchException as e:
            revcom_match += 1
            halo_index, molecular_tag = e[0], e[1]
        except NoIndexMatchException as e:
            non_match += 1
            halo_index, molecular_tag = e[0], e[1]
        # TODO modify headers of the reads and split by index either in dict, yield, or write to file
        # possibly use scilifelab.utils.fastq_utils.FastQWriter()
        if halo_index and molecular_tag:
            for read_num, read in enumerate([read_1, read_2]):
                read[0] = read[0] + "{}:{}:".format(halo_index, molecular_tag)
                file_name = os.path.join(output_dir, "{}_R{}.fastq".format(halo_index, read_num+1))
                with open(file_name, 'a') as f:
                    f.write("\n".join(read))
    import pdb; pdb.set_trace()


def parse_index(index_seq, halo_index_dict=None, halo_index_revcom_dict=None, halo_index_length=None, molecular_tag_length=None):
    """
    Split an index up into its haloplex index and the random molecular tag.
    Returns the halo index and the molecular tag as a tuple.
    """
    halo_index, molecular_tag = index_seq[:halo_index_length], index_seq[halo_index_length:molecular_tag_length]

    # TODO add check against known indexes (Bio.pairwise2) -- this will also determine default molecular tag length
    #from Bio import pairwise2 as pw2
    #    pw2.align.globalms(halo_index, key, 2, -1, -1, -1)
    if not halo_index in halo_index_dict.keys():
        if halo_index in halo_index_revcom_dict.keys():
            raise RevComMatchException(None, halo_index, molecular_tag)
        else:
            raise NoIndexMatchException(None, halo_index, molecular_tag)
    #    return None, None, None
    #index_name = halo_index_dict[halo_index]
    #return halo_index, molecular_tag, index_name
    return halo_index, molecular_tag


def load_indexes(csv_file):
    """
    Load known indexes from a csv file.
    csv file should be in format:
        sequence[,index_name]
    where index_name is optional.
    Returns a dict of sequence->name pairs.
    """
    index_dict          = {}
    index_dict_revcom   = {}
    with open(csv_file, 'r') as f:
        for line in f:
            # could also use csv.sniffer to dynamically determine delimiter
            index = re.split(r'[\t,;]', line.strip())
            # include reverse complement
            rev_com_index = Bio.Seq.Seq(index[0]).reverse_complement().tostring()
            try:
                index_dict[index[0]]                = index[1]
                index_dict_revcom[rev_com_index]    = index[1]
            except IndexError:
                index_dict[index[0]]                = None
                index_dict_revcom[rev_com_index]    = None
    return index_dict, index_dict_revcom

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample-sheet", help="Sample sheet describing the indexes used during library preparation.")
    parser.add_argument("-i", "--halo-index-file", required=True, help="File containing haloplex indexes (one per line, optional name in second column separated by tab, comma, or semicolon).")
    # TODO make this not required by dynamically matching to indexes
    parser.add_argument("-l", "--halo-index-length", type=int, required=True, help="The length of the haloplex index. Required if indexes are not supplied (option -i).")
    parser.add_argument("-m", "--molecular-tag-length", type=int, help="The length of the (random) molecular tag. If not specified, the remainder of the read after the halo index is used.")
    #parser.add_argument("-r", "--read-index", dest="read_index", type=int, help="Which read is the index.")
    #parser.add_argument("-s", "--single-read", dest="single_read", action="store_true", help="Specify that the data is single-read (not paired-end); default false.")
    #parser.add_argument("-d", "--directory", help="Directory containing demultiplexed data.")
    parser.add_argument("-1", "--read-one", required=True, help="Read 1 fastq file.")
    parser.add_argument("-2", "--read-two", required=True, help="Read 2 fastq file.")
    parser.add_argument("-r", "--read-index", required=True, help="Index read fastq file.")

    arg_vars = vars(parser.parse_args())

    main(**arg_vars)
