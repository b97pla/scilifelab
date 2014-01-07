"""
Demultiplex haloplex data including molecular tags.
"""
# TODO add pairwise alignment of indexes to correct for sequencing error (biopython's Bio.pairwise2)

from __future__ import print_function

import argparse
#import collections
import itertools
import os
import re
import shutil
import subprocess
import sys

from Bio import Seq, pairwise2
#from scilifelab.illumina.hiseq import HiSeqSampleSheet
#from scilifelab.utils.fastq_utils import FastQParser

# TODO add index file reading
# TODO switch from open() to subprocess pipes?
# TODO add directory processing
# TODO add status updates (every 10K read or whatever)
# TODO add minimum read cutoffs (i.e. drop indexes with reads fewer than X)


def main(read_one, read_two, read_index, data_directory, read_index_num, output_directory, halo_index_file, halo_index_length, molecular_tag_length=None, force_overwrite=False):

    if not output_directory:
        raise SyntaxError("Must specify output directory.")
    else:
        output_directory = create_output_dir(output_directory, force_overwrite)

    if halo_index_file:
        # Load indexes, names from a file
        index_dict, index_revcom_dict = load_index_file(halo_index_file)
    elif halo_index_length:
        assert(type(halo_index_length) == int and halo_index_length > 0), "Haloplex index length must be a positive integer."
        if molecular_tag_length:
            assert(type(molecular_tag_length) == int and molecular_tag_length >= 0), "Molecular tag length must be a postive integer."
        # else indexes will be extracted from the data itself based on length
        index_dict, index_revcom_dict = None, None
    else:
        raise SyntaxError("Either an index file or the index length must be specified.")

    if read_one and read_two and read_index:
        if data_directory or read_index_num:
            raise SyntaxError("Ambiguous: too many options specified. Specify either file paths or directory and read index number.")
        else:
            if index_dict:
                # Indexes are user-supplied
                reads_processed, num_match, num_nonmatch = parse_readset_byindexdict(read_one, read_two, read_index, index_dict, index_revcom_dict, output_directory)
                print("\nInfo: Processing complete; {} reads processed ({} matching, {} non-matching).".format(reads_processed, num_match, num_nonmatch), file=sys.stderr)
            else:
                # Indexes will be pulled from data itself based on user-supplied length
                reads_processed = parse_readset(read_one, read_two, read_index, output_directory, halo_index_length, molecular_tag_length)
                print("\nInfo: Processing complete; {} reads processed..".format(reads_processed), file=sys.stderr)
    elif data_directory and read_index_num:
        parse_directory(data_directory, read_index_num)
        #for readset in parse_directory(directory, read_index_num):
        # etc.
        # requires changing parse_directory to generator and changing parse_readset to accept a dict for the reads e.g. {1:x, 2:x, i:x}
        #parse_directory(directory, read_index_num, index_dict, index_revcom_dict, output_directory, halo_index_length, molecular_tag_length)
    else:
        raise SyntaxError("Either a directory and read index number or explicit paths to sequencing files must be specified.")



# TODO turn this into a generator and use it in a for loop above
def parse_directory(data_directory, read_index_num):
    """
    Searches the directory for fastq file sets and calls parse_readset() on them.
    """
    raise NotImplementedError("I haven't implemented this yet, so don't go using it.")

# TODO parallelize
# TODO add checking against reverse_complement sequences
# TODO there's probably a good way to combine these two parse_readset functions
def parse_readset_byindexdict(read_1_fq, read_2_fq, read_index_fq, index_dict, index_revcom_dict, output_directory, verbose=True):
    """
    Parse input fastq files, searching for matches to each index.
    """

    fqp_1, fqp_2, fqp_ind = map(FastQParser, (read_1_fq, read_2_fq, read_index_fq))

    #match, non_match, revcom_match = 0, 0, 0
    reads_processed, num_match, num_nonmatch = 0, 0, 0

    if verbose:
        total_lines_in_file = sum(1 for line in open(read_1_fq))

    for read_1, read_2, read_ind in itertools.izip(fqp_1, fqp_2, fqp_ind):
        read_ind_seq = read_ind[1]
        for supplied_index in index_dict.keys():
            if re.match(supplied_index, read_ind_seq):
                index_len = len(supplied_index)
                index, molecular_tag = read_ind_seq[:index_len], read_ind_seq[index_len:]
                num_match += 1

                read_list = [read_1, read_2]
                modify_reads(read_list, index, molecular_tag)
                sample_name = index_dict[supplied_index] if index_dict[supplied_index] else supplied_index
                write_reads_to_disk(read_list, sample_name, output_directory)

                break
            else:
                num_nonmatch += 1
        reads_processed += 1

        if total_lines_in_file:
            percentage = 100 * ((reads_processed * 4.0) / total_lines_in_file)
            if percentage % 1 == 0:
                print_progress(percentage)

    print("\nProcessed {} reads, processing complete.".format(reads_processed), file=sys.stderr)
    if num_match == 0:
        print("\nWarning: no reads matched supplied indexes; nothing written.", file=sys.stderr)
    return reads_processed, num_match, num_nonmatch

def print_progress(percentage_complete):
    """
    Prints a little progress bar on the current line after clearing it.
    """
    _, term_cols = os.popen('stty size', 'r').read().split()
    term_cols = int(term_cols) - 10
    progress = int((term_cols * percentage_complete / 100))
    sys.stderr.write('\r')
    sys.stderr.write("[{progress:<{cols}}] {percentage}%".format(progress="="*progress, cols=term_cols, percentage=percentage_complete))
    sys.stderr.flush()


# TODO there's probably a good way to combine these two parse_readset functions
def parse_readset(read_1_fq, read_2_fq, read_index_fq, output_directory, halo_index_length, molecular_tag_length=None, verbose=True):
    """
    Parse input fastq files by index reads.
    """
    fqp_1, fqp_2, fqp_ind = map(FastQParser, (read_1_fq, read_2_fq, read_index_fq))

    #match, non_match, revcom_match = 0, 0, 0
    reads_processed, reads_match, reads_nonmatch = 0, 0, 0

    if verbose:
        total_lines_in_file = sum(1 for line in open(read_1_fq))

    for read_1, read_2, read_ind in itertools.izip(fqp_1, fqp_2, fqp_ind):
        read_ind_seq = read_ind[1]
        index, molecular_tag = read_ind_seq[:halo_index_length], read_ind_seq[halo_index_length:molecular_tag_length] 

        read_list = [read_1, read_2]
        modify_reads(read_list, index, molecular_tag)
        sample_name = index
        write_reads_to_disk(read_list, sample_name, output_directory)

        reads_processed += 1

        if total_lines_in_file:
            percentage = 100 * ((reads_processed * 4.0) / total_lines_in_file)
            if percentage % 1 == 0:
                print_progress(percentage)

    return reads_processed

def modify_reads(read_list, index, molecular_tag):
    """
    Add the extra tags to the header of each read
    """
    for read in read_list:
        read[0] = read[0] + "{}:{}:".format(index, molecular_tag)

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

# TODO see if it's faster to keep a dict of { index -> filehandle }
def write_reads_to_disk(read_list, sample_name, output_directory):
    """
    Write a Fastq read to the appropriate file.
    """
    for read_num, read in enumerate(read_list):
        file_name = "{sample_name}_R{read_num}.fastq".format(sample_name=sample_name, read_num=read_num+1)
        file_name = os.path.join(output_directory, file_name)
        with open(file_name, 'a') as f:
            f.write("\n".join(read) + "\n")


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
            rev_com_index = Seq.Seq(index[0]).reverse_complement().tostring()
            try:
                index_dict[index[0]]                = index[1]
                index_dict_revcom[rev_com_index]    = index[1]
            except IndexError:
                index_dict[index[0]]                = None
                index_dict_revcom[rev_com_index]    = None
    return index_dict, index_dict_revcom

# I'm just putting this class here temporarily so I can test this without needing to load in all the other scilifelab jazz
class FastQParser(object):
    """Parser for fastq files, possibly compressed with gzip.
       Iterates over one record at a time. A record consists
       of a list with 4 elements corresponding to 1) Header,
       2) Nucleotide sequence, 3) Optional header, 4) Qualities"""

    def __init__(self, file, filter=None):
        self.fname = file
        self.filter = filter

        if file.endswith(".gz"):
            self._fh = subprocess.Popen(["gunzip", "-d", "-c", file], stdout = subprocess.PIPE, bufsize = 1).stdout
        else:
            self._fh = subprocess.Popen(["cat", file], stdout = subprocess.PIPE, bufsize = 1).stdout
        self._records_read = 0
        self._next = self.setup_next()

    def __iter__(self):
        return self

    def next(self):
        return self._next(self)

    def setup_next(self):
        """Return the function to return the next record
        """
        if self.filter is None or len(self.filter.keys()) == 0:
            def _next(self):
                self._records_read += 1
                return [self._fh.next().strip() for n in range(4)]
        else:
            def _next(self):
                while True:
                    record = [self._fh.next().strip() for n in range(4)]
                    header = parse_header(record[0])
                    skip = False
                    for k, v in self.filter.items():
                        if k in header and header[k] not in v:
                            skip = True
                            break
                    if not skip:
                        self._records_read += 1
                        return record
        return _next

    def name(self):
        return self.fname

    def rread(self):
        return self._records_read

    def seek(self,offset,whence=None):
        self._fh.seek(offset,whence)

    def close(self):
        self._fh.close()

def parse_header(header):
    """Parses the FASTQ header as specified by CASAVA 1.8.2 and returns the fields in a dictionary
       @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
    """
    if header[0] != '@':
        return None

    instrument, run_number, flowcell_id, lane, tile, x_pos, y_pos_read, is_filtered, control_number, index = header[1:].split(":")
    y_pos, read = y_pos_read.split()
    return {'instrument': str(instrument),
            'run_number': int(run_number),
            'flowcell_id': str(flowcell_id),
            'lane': int(lane),
            'tile': int(tile),
            'x_pos': int(x_pos),
            'y_pos': int(y_pos),
            'read': int(read),
            'is_filtered': (is_filtered == 'Y'),
            'control_number': int(control_number),
            'index': str(index)} # Note that MiSeq Reporter outputs a SampleSheet index rather than the index sequence

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output-directory", required=True,
                                help="The directory to be used for storing output data. Required.")

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
                                help="Directory containing fastq read data. Requires specifying read number (-n).")
    parser.add_argument("-n", "--read-index-num", type=int,
                                help="Which read is the index (e.g. 1, 2, 3).")

    parser.add_argument("-f", "--force-overwrite", action="store_true",
                                help="Force overwrite to output directory.")

    # TODO parser.add_argument("-s", "--single-read", action="store_true", help="Specify that the data is single-read (not paired-end). Default false.")

    arg_vars = vars(parser.parse_args())

    main(**arg_vars)
