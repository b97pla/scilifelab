"""
Demultiplex haloplex data including molecular tags.
"""
# TODO add pairwise alignment of indexes to correct for sequencing error (biopython's Bio.pairwise2)

from __future__ import print_function

import argparse
import Bio
#import collections
import itertools
import os
import re
import shutil
import sys

#from scilifelab.illumina.hiseq import HiSeqSampleSheet
from scilifelab.utils.fastq_utils import FastQParser


# TODO these classes don't print the correct error message - maybe a problem with super? review docs
class NoIndexMatchException(Exception):
    def __init__(self, message, *args):
        if not message:
            message = "no match to supplied indexes."
        super(NoIndexMatchException, self).__init__(self, message)
        #Exception.__init__(self, message)
        self.args = args

class RevComMatchException(Exception):
    def __init__(self, message, *args):
        if not message:
            message = "no match to supplied indexes but match to reverse complement of supplied indexes. Double-check read direction."
        super(RevComMatchException, self).__init__(self, message)
        Exception.__init__(self, message)
        self.args = args



def main(read_one, read_two, read_index, data_directory, read_index_num, output_directory, halo_index_file, halo_index_length, molecular_tag_length=None, force_overwrite=False):
#def main(**kwargs):

    if not output_directory:
        raise SyntaxError("Must specify output directory.")
    output_directory = create_output_dir(output_directory, force_overwrite)

    # Check for minimum required index information from user
    if halo_index_file:
        index_dict, index_revcom_dict = load_index_file(halo_index_file)
    elif halo_index_length:
        assert(type(halo_index_length) == int and halo_index_length > 0), "Haloplex index length must be a positive integer."
        if molecular_tag_length:
            assert(type(molecular_tag_length) == int and molecular_tag_length >= 0), "Molecular tag length must be a postive integer."
        index_dict, index_revcom_dict = None, None
    else:
        raise SyntaxError("Either an index file or the index length must be specified.")

    # Check for minimum required fastq location data from user and process
    if read_one and read_two and read_index:
        if data_directory or read_index_num:
            raise SyntaxError("Ambiguous: too many options specified. Specify either file paths or directory and read index number.")
        else:
            parse_readset(read_one, read_two, read_index, index_dict, index_revcom_dict, output_directory, halo_index_length, molecular_tag_length)
    elif data_directory and read_index_num:
        parse_directory(data_directory, read_index_num)
        #for readset in parse_directory(directory, read_index_num):
        # etc.
        # requires changing parse_directory to generator and changing parse_readset to accept a dict for the reads e.g. {1:x, 2:x, i:x}
        #parse_directory(directory, read_index_num, index_dict, index_revcom_dict, output_directory, halo_index_length, molecular_tag_length)
    else:
        raise SyntaxError("Either a directory and read index number or explicit paths to sequencing files must be specified.")
        print("End of function", file=sys.stderr)


def create_output_dir(output_directory, force_overwrite):
    """
    Create the output directory, forcing overwrite if the -f flag is passed; otherwise, fail if file exists."
    """
    output_directory = os.path.abspath(output_directory)
    if os.path.exists(output_directory):
        if force_overwrite:
            print("Warning: removing/overwriting output directory \"{}\"...".format(output_directory), file=sys.stderr, end="")
            shutil.rmtree(output_directory)
            print(" removed.", file=sys.stderr)
    # This can be outside the conditional because I don't want to raise my own exception in case of failure
    os.makedirs(output_directory)
    return output_directory

# TODO turn this into a generator and use it in a for loop above
#def parse_directory(data_directory, read_index_num, index_dict, index_revcom_dict, halo_index_length, molecular_tag_length):
def parse_directory(data_directory, read_index_num):
    """
    Searches the directory for fastq file sets and calls parse_readset() on them.
    """
    raise NotImplementedError("I haven't implemented this yet, so don't go using it.")

def parse_readset(read_1_fq, read_2_fq, read_index_fq, index_dict, index_revcom_dict, output_directory, halo_index_length, molecular_tag_length=None, force_append=False):
    """
    Parse input fastq files by index reads.
    """

    fqp_1, fqp_2, fqp_ind = map(FastQParser, (read_1_fq, read_2_fq, read_index_fq))

    match, non_match, revcom_match = 0, 0, 0

    for read_1, read_2, read_ind in itertools.izip(fqp_1, fqp_2, fqp_ind):
        # Get the haloplex index and the molecular tag from the index read's sequence data
        try:
            halo_index, molecular_tag = parse_index(read_ind[1], index_dict, index_revcom_dict,
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
                file_name = os.path.join(output_directory, "{}_R{}.fastq".format(halo_index, read_num+1))
                with open(file_name, 'a') as f:
                    f.write("\n".join(read) + "\n")


def parse_index(index_seq, index_dict=None, index_revcom_dict=None, halo_index_length=None, molecular_tag_length=None):
    """
    Split an index up into its haloplex index and the random molecular tag.
    Returns the halo index and the molecular tag as a tuple.
    """
    halo_index, molecular_tag = index_seq[:halo_index_length], index_seq[halo_index_length:molecular_tag_length]

    # TODO add check against known indexes (Bio.pairwise2) -- this will also determine default molecular tag length
    #from Bio import pairwise2 as pw2
    #    pw2.align.globalms(halo_index, key, 2, -1, -1, -1)
    #if not halo_index in index_dict.keys():
    #    if halo_index in index_revcom_dict.keys():
    #        raise RevComMatchException(None, halo_index, molecular_tag)
    #    else:
    #        raise NoIndexMatchException(None, halo_index, molecular_tag)
    #    return None, None, None
    #index_name = index_dict[halo_index]
    #return halo_index, molecular_tag, index_name
    return halo_index, molecular_tag


def load_index_file(csv_file):
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

def check_input(**kwa):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output-directory", required=True,
                                help="The directory to be used for storing output data.")

    parser.add_argument( "-i", "--halo-index-file",
                                help="File containing haloplex indexes (one per line, optional name " \
                                     "in second column separated by tab, comma, or semicolon).")
    parser.add_argument( "-l", "--halo-index-length", type=int,
                                help="The length of the haloplex index. Required if indexes are not supplied (option -i).")
    parser.add_argument( "-m", "--molecular-tag-length", type=int,
                                help="The length of the (random) molecular tag. If not specified, " \
                                     "the remainder of the read after the halo index is used.")

    parser.add_argument("-1", "--read-one",
                                help="Read 1 fastq file.")
    parser.add_argument("-2", "--read-two",
                                help="Read 2 fastq file.")
    parser.add_argument("-r", "--read-index",
                                help="Index read fastq file.")

    parser.add_argument("-d", "--data-directory",
                                help="Directory containing fastq read data.")
    parser.add_argument("-n", "--read-index-num", type=int,
                                help="Which read is the index (e.g. 1, 2, 3).")

    parser.add_argument("-f", "--force-overwrite", action="store_true",
                                help="Force overwrite to output directory.")

    # TODO parser.add_argument("-s", "--single-read", action="store_true", help="Specify that the data is single-read (not paired-end). Default false.")

    # TODO ensure exclusive subcommands for indexes, fastq data: check namespaces
    # TODO set unset variables to None? or check for presence of variables here
    arg_vars = vars(parser.parse_args())

    main(**arg_vars)
