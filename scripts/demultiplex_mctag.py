"""
Demultiplex haloplex data including molecular tags.
"""

from __future__ import print_function

import argparse

from scilifelab.utils.fastq_utils import FastQParser


def main(read_one, read_two, read_ind, sample_sheet, halo_index_length, molecular_tag_length):

    parse_readset(read_one, read_two, read_ind)


def parse_readset(read_one, read_two, read_ind):
    """
    Parse input fastq files by index reads.
    """
    fp1, fp2, fpind = map(FastQParser, (read_one, read_two, read_ind))
    for r1, r2, rind in fp1, fp2, fpind:
        print("\n".join(r1, r2, rind))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--sample-sheet", help="Sample sheet describing the indexes used during library preparation.")
    parser.add_argument("-l", "--halo-index-length", help="The length of the haloplex index.")
    parser.add_argument("-m", "--molecular-tag-length", help="The length of the (random) molecular tag.")
    #parser.add_argument("-r", "--read-index", dest="read_index", type=int, help="Which read is the index.")
    #parser.add_argument("-s", "--single-read", dest="single_read", action="store_true", help="Specify that the data is single-read (not paired-end); default false.")
    #parser.add_argument("-d", "--directory", help="Directory containing demultiplexed data.")
    parser.add_argument("-1", "--read-one", help="Data read 1")
    parser.add_argument("-2", "--read-two", help="Data read 2")
    parser.add_argument("-i", "--read-index", help="Index read.")

    arg_vars = vars(parser.parse_args())

    main(**arg_vars)
