"""
Demultiplex haloplex data including molecular tags.
"""

from __future__ import print_function

import argparse
import collections
import datetime
import fcntl
import itertools
import os
import random
import re
import resource
import shlex
import shutil
import subprocess
import sys
import time

#from Bio import Seq, pairwise2
#from scilifelab.utils.fastq_utils import FastQParser, FastQWriter

# TODO memoize sequence corrections? optimize somehow if possible
# TODO ensure read 1,2 files are paired (SciLifeLab code)
# TODO add directory processing


def main(read_one, read_two, read_index, data_directory, read_index_num, output_directory, index_file, max_mismatches=1, force_overwrite=False, progress_interval=1000):
    check_input(read_one, read_two, read_index, data_directory, read_index_num, output_directory, index_file, max_mismatches, progress_interval)
    output_directory = create_output_dir(output_directory, force_overwrite)
    # TODO NOTE: I changed this just for Jimmie's run
    index_dict = load_index_file(index_file, usecols=(2,1))
    check_index_distances(index_dict.keys(), max_mismatches)
    start_time = datetime.datetime.now()
    if read_one and read_two and read_index:
        reads_processed, num_match, num_ambigmatch, num_nonmatch, num_corrected = parse_readset_byindexdict(read_one, read_two, read_index, index_dict, output_directory, max_mismatches, progress_interval)
    else:
        #for readset in parse_directory(data_directory, read_index_num):
            #read_one, read_two, read_index = readset
        reads_processed, num_match, num_ambigmatch, num_nonmatch, num_corrected = parse_readset_byindexdict(read_one, read_two, read_index, index_dict, output_directory)
    elapsed_time = time.strftime('%H:%M:%S', time.gmtime((datetime.datetime.now() - start_time).total_seconds()))
    print(  "\nProcessing complete in {elapsed_time}:\n\t" \
            "{reads_processed} reads processed\n\t" \
            "{num_match:>{pad_length}} ({num_match_percent:>6.2f}%) matched to supplied indexes\n\t" \
            "{num_corrected:>{pad_length}} ({num_corrected_percent:>6.2f}%) corrected indexes\n\t" \
            "{num_ambigmatch:>{pad_length}} ({num_ambigmatch_percent:>6.2f}%) matches to more than one supplied index.\n\t" \
            "{num_nonmatch:>{pad_length}} ({num_nonmatch_percent:>6.2f}%) unmatched to supplied indexes".format(
                elapsed_time=elapsed_time,
                reads_processed=reads_processed, num_match=num_match, num_nonmatch=num_nonmatch, num_ambigmatch=num_ambigmatch, num_corrected=num_corrected,
                num_match_percent       = (100.0 * num_match)/reads_processed,
                num_corrected_percent   = (100.0 * num_corrected)/reads_processed,
                num_ambigmatch_percent  = (100.0 * num_ambigmatch)/reads_processed,
                num_nonmatch_percent    = (100.0 * num_nonmatch)/reads_processed,
                pad_length = len(str(reads_processed))), file=sys.stdout)


def check_input(read_one, read_two, read_index, data_directory, read_index_num, output_directory, index_file, max_mismatches, progress_interval):
    """
    Check user-supplied inputs for validity, completeness.
    """
    if not output_directory:
        raise SyntaxError("Must specify output directory.")
    if not index_file:
        raise SyntaxError("Must specify file containing indexes.")
    if not (read_one and read_index):
        raise SyntaxError("Must speify both data and index reads.")
    if (read_one or read_two or read_index) and (data_directory or read_index_num):
        raise SyntaxError("Ambiguous: too many options specified. Specify either file paths or directory and read index number.")
    if not (read_one and read_two and read_index) or (data_directory and read_index_num):
        raise SyntaxError("Insufficient information: either a directory and read index number or explicit paths to sequencing files must be specified.")
    try:
        assert(type(progress_interval) == int and progress_interval > 0)
    except AssertionError:
        raise SyntaxError("Progress interval must be a positive integer.")
    try:
        assert(type(max_mismatches) == int and max_mismatches >= 0)
    except AssertionError:
        raise SyntaxError("Maximum mismatches in error correction must be >= 0.")


def parse_directory(data_directory, read_index_num):
    """
    Searches the directory for fastq file sets and calls parse_readset() on them.
    """
    raise NotImplementedError("I haven't implemented this yet, so don't go using it.")
    # possibly implement as generator, calling parse_readset_byindexdict in a for loop from the calling loop


def parse_readset_byindexdict(read_1_fq, read_2_fq, read_index_fq, index_dict, output_directory, max_mismatches=1, progress_interval=1000):
    """
    Parse input fastq files, searching for matches to each index.
    """
    print("Processing read set associated with \"{}\" using user-supplied indexes.".format(read_1_fq), file=sys.stderr)
    print("Maximum number of mismatches for error correction is {}.".format(max_mismatches), file=sys.stderr)
    reads_processed, num_match, num_ambigmatch, num_nonmatch, num_corrected = 0, 0, 0, 0, 0
    fqp_1, fqp_2, fqp_ind = map(FastQParser, (read_1_fq, read_2_fq, read_index_fq))
    print("Counting total number of lines in fastq files...", file=sys.stderr, end="")
    # I think du -k * 16 / 1.024 should give approximately the right number for any number of reads greater than 1000 or so
    total_lines_in_file = int(subprocess.check_output(shlex.split("wc -l {}".format(read_1_fq))).split()[0])
    print(" complete.", file=sys.stderr)
    if not progress_interval: progress_interval = 1000
    if progress_interval > (total_lines_in_file / 4): progress_interval = (total_lines_in_file / 4)
    index_fh_dict = collections.defaultdict(list)
    print("Demultiplexing...", file=sys.stderr)
    time_started = datetime.datetime.now()
    for read_1, read_2, read_ind in itertools.izip(fqp_1, fqp_2, fqp_ind):
        read_ind_seq = read_ind[1]
        matches_dict = collections.defaultdict(list)
        # Sort indexes by descending length to match longer indexes first
        for supplied_index in sorted(index_dict.keys(), key=lambda x: (-len(x))):
            mismatches = find_dist(supplied_index, read_ind_seq, max_mismatches)
            matches_dict[mismatches].append(supplied_index)
            if mismatches == 0:
                break
        for x in range(0, max_mismatches+1):
            if matches_dict.get(x):
                if len(matches_dict.get(x)) == 1:
                    # Single unamibiguous match
                    index_seq       = matches_dict[x][0]
                    index_len       = len(index_seq)
                    molecular_tag   = read_ind_seq[index_len:]
                    modify_reads( (read_1, read_2), index_seq, molecular_tag)
                    sample_name     = index_dict[index_seq] if index_dict[index_seq] else index_seq
                    data_write_loop(read_1, read_2, sample_name, output_directory, index_fh_dict, index_seq)
                    num_match      += 1
                    if not x == 0:
                        num_corrected += 1
                    break
                else:
                    # Ambiguous match
                    sample_name     = "Ambiguous"
                    index_seq_list  = ",".join(matches_dict.get(x))
                    modify_reads( (read_1, read_2), index_seq_list, read_ind_seq)
                    data_write_loop(read_1, read_2, sample_name, output_directory, index_fh_dict, sample_name)
                    num_ambigmatch += 1
                    break
        else:
            # No match
            sample_name     = "Undetermined"
            modify_reads( (read_1, read_2), "", read_ind_seq)
            data_write_loop(read_1, read_2, sample_name, output_directory, index_fh_dict, sample_name)
            num_nonmatch   += 1
        reads_processed    += 1
        if reads_processed % progress_interval == 0:
            print_progress(reads_processed, (total_lines_in_file / 4), time_started=time_started)
    return reads_processed, num_match, num_ambigmatch, num_nonmatch, num_corrected


def data_write_loop(read_1, read_2, sample_name, output_directory, index_fh_dict, index):
    """
    Writes data using FastQParser, closing files if we open too many.
    """
    for read_num, read in enumerate([read_1, read_2]):
        try:
            index_fh_dict[index][read_num].write(read)
        except IndexError:
            file_path   = os.path.join(output_directory, "{sample_name}_R{read_num}.fastq".format(sample_name=sample_name, read_num=read_num+1))
            try:
                index_fh_dict[index].append(FastQWriter(file_path))
            except IOError as e:
                # Too many open filehandles
                if e.errno == 24:
                    #print("Warning: too many open file handles. Closing...", file-sys.stderr)
                    for fh1, fh2 in index_fh_dict.values():
                        map(file.close, [ fh1, fh2 ] )
                    index_fh_dict[index].append(FastQWriter(file_path))
                    index_fh_dict[index][read_num].write(read)
                else:
                    raise IOError(e)
            except ValueError:
                # File was closed previously
                index_fh_dict[index][read_num].reopen()

# TODO compare to Bio.align.pairwise2 for speed
# TODO possibly @memoize somehow
def find_dist(str_01, str_02, max_mismatches=None, approach="shorten"):
    """
    Find the number of mismatches between two strings. The longer string is truncated to the length of the shorter unless approach "lengthen".
    """
    if len(str_01) > len(str_02):
        if approach == "lengthen":
            str_02 = "{:<{length}}".format(str_02, length=len(str_01))
        else:
            str_01 = str_01[:len(str_02)]
    elif len(str_02) > len(str_01):
        if approach == "lengthen":
            str_01 = "{:<{length}}".format(str_01, length=len(str_02))
        else:
            str_02 = str_02[:len(str_01)]
    mismatches = 0
    for a, b in itertools.izip(str_01, str_02):
        if a != b:
            mismatches += 1
            if max_mismatches and mismatches > max_mismatches:
                break
    return mismatches


def print_progress(processed, total, type='text', time_started=None, leading_text=""):
    """
    Prints the progress, either in text or in visual form.
    """
    percentage_complete = float(processed) / total
    if time_started:
        completion_time = estimate_completion_time(time_started, percentage_complete)
    else:
        time_started = "-"
    _, term_cols = os.popen('stty size', 'r').read().split()
    term_cols = int(term_cols) - 10 - len(leading_text)
    sys.stderr.write('\r')
    if type == 'bar' and term_cols > 10:
        progress = int(term_cols * percentage_complete)
        sys.stderr.write("{leading_text}[{progress:<{cols}}] {percentage:0.2f}%".format(
            leading_text=leading_text, progress="="*progress, cols=term_cols, percentage=percentage_complete*100))
    else:
        sys.stderr.write("{leading_text}{processed}/{total} items processed ({percentage_complete:0.2f}% finished) ETA: {completion_time}".format(
            leading_text=leading_text, processed=processed, total=total, percentage_complete=percentage_complete*100, completion_time=completion_time))
    sys.stderr.flush()

def estimate_completion_time(start_time, percent_complete):
    """
    http://xkcd.com/612/
    """
    if not type(start_time) == datetime.datetime:
        return None
    seconds_elapsed = (datetime.datetime.now() - start_time).total_seconds()
    seconds_total = seconds_elapsed / percent_complete
    seconds_left = seconds_total - seconds_elapsed
    # More than a day remaining
    if seconds_left > 86400:
        days    = int(seconds_left // 86400)
        seconds = seconds_left - (days * 86400)
        return "{}:{}".format(days, time.strftime('%H:%M:%S', time.gmtime(seconds)))
    else:
        return time.strftime('%H:%M:%S', time.gmtime(seconds_left))


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
    os.makedirs(output_directory)
    return output_directory


def load_index_file(csv_file, usecols=(0,1)):
    """
    Load known indexes from a csv file.
    csv file should be in format:
        sequence[,index_name]
    where index_name is optional.
    Returns a dict of sequence->name pairs.
    """
    if not len(usecols) == 2:
        raise ValueError("Warning: only two columns (index, sample_name) can be selected.", file=sys.stderr)
    if not usecols == (0,1):
        print("Warning: using non-standard columns for index, sample_name ({} instead of (0,1))".format(usecols), file=sys.stderr)
    index_column, name_column = usecols
    index_dict          = {}
    with open(csv_file, 'r') as f:
        for line in f:
            # could also use csv.sniffer to dynamically determine delimiter
            index_line = re.split(r'[\t,;]', line.strip())
            try:
                index_dict[ index_line[index_column] ] = index_line[name_column]
            except IndexError:
                index_dict[ index_line[index_column] ] = None
    return index_dict


def check_index_distances(index_list, max_mismatches):
    """
    Determines if too many mismatches are allowed for this set of indexes to resolve unambiguously.
    """
    for i1, i2 in  itertools.combinations(index_list, r=2):
        if find_dist(i1, i2, max_mismatches, approach="lengthen") <= max_mismatches:
            print("Warning: indexes \"{}\" and \"{}\" are insufficiently different for the specified number of mismatches ({}). Reads matching either index will be classified as ambiguous.".format(i1, i2, max_mismatches), file=sys.stderr)


###

# TODO This doesn't really belong here and should probably be its own module
def count_top_indexes(count_num, index_file, index_length, progress_interval):
    """
    Determine the most common indexes, sampling at most 200,000 reads.
    """
    assert(type(count_num) == int and count_num > 0), "Number passed must be a positive integer."
    fqp_ind = FastQParser(index_file)
    # This should perhaps be added to the FastQParser class
    print("Counting total number of lines in fastq file...", file=sys.stderr, end="")
    total_lines = int(subprocess.check_output(shlex.split("wc -l {}".format(index_file))).split()[0])
    total_reads = total_lines / 4
    print(" complete.", file=sys.stderr)
    index_tally = collections.defaultdict(int)
    reads_processed = 0
    # Subsample if file is large
    if (total_reads) > 200000:
        print("Subsampling 200,000 reads from index file...", file=sys.stderr)
        fqp_ind = iter_sample_fast(fqp_ind, 200000, total_reads)
        print("Complete.", file=sys.stderr)
        total_reads = 200000
    print("Tallying indexes in {} records...".format(total_reads), file=sys.stderr)
    start_time = datetime.datetime.now()
    for index in fqp_ind:
        index_read_seq = index[1]
        index_seq = index_read_seq[:index_length]
        index_tally[index_seq] += 1
        reads_processeds += 1
        if reads_processed % progress_interval == 0:
            print_progress(reads_processed, total_reads, start_time)
    print("\n", file=sys.stderr)
    if count_num > len(index_tally.keys()):
        print("Number of indexes found ({}) is fewer than those requested ({}). Printing all indexes found.".format(len(index_tally.keys()), count_num), file=sys.stderr)
        print("Printing indexes...", file=sys.stderr())
        count_num = len(index_tally.keys())
    print("{:<20} {:>20} {:>11}".format("Index", "Occurences", "Percentage"))
    for index, _ in sorted(index_tally.items(), key=(lambda x: x[1]), reverse=True)[:count_num]:
        percentage = (100.0 * index_tally[index] ) / total_reads
        print("{:<20} {:>20,} {:>10.2f}%".format(index, index_tally[index], percentage))

def iter_sample_fast(iterable, samplesize, total_size):
    """
    http://stackoverflow.com/questions/12581437/python-random-sample-with-a-generator/12583436#12583436
    """
    results = []
    iterator = iter(iterable)
    # Fill in the first samplesize elements:
    try:
        for _ in xrange(samplesize):
            results.append(iterator.next())
            print_progress(len(results), 200000)
    except StopIteration:
        raise ValueError("Sample larger than population.")
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items
            print_progress(i, total_size)
    return results



# SCILIFELAB CODE
# I'm just putting this class here temporarily so I can test this without needing to load in all the other scilifelab jazz
# TODO I wonder if I can adjust this to use a larger byte size instead of scanning the lines as it probably does now -- the increased read block size would help speed enormously
class FastQParser(object):
    """
    Parser for fastq files, possibly compressed with gzip.
    Iterates over one record at a time. A record consists
    of a list with 4 elements corresponding to 1) Header,
    2) Nucleotide sequence, 3) Optional header, 4) Qualities
    """

    def __init__(self, file, filter=None):
        self.fname = file
        self.filter = filter

        if file.endswith(".gz"):
            self._fh = subprocess.Popen(["gunzip", "-d", "-c", file], stdout = subprocess.PIPE, bufsize = 1).stdout
        else:
            # TODO I'm not sure this is faster than just open(file, 'w')
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

class FastQWriter:

    def __init__(self,file):
        self.fname = file
        if file.endswith(".gz") or file.endswith(".gzip"):
            command = "gzip -c  >> {}".format(file)
            self._fh = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE).stdin
        else:
            self._fh = open(file, 'a')
            #self._fh = subprocess.Popen("cat >> {}".format(file), shell=True, stdin=subprocess.PIPE).stdin
            #self._fh = subprocess.Popen(open(file, 'a+'), stdin=subprocess.PIPE)
        self._records_written = 0

    def name(self):
        return self.fname

    def write(self,record):
        self._fh.write("{}\n".format("\n".join([r.strip() for r in record])))
        self._records_written += 1

    def rwritten(self):
        return self._records_written

    def close(self):
        #self.stdin.flush()
        self._fh.close()

    def reopen(self):
        _records_written = self._records_written
        self.__init__(self.fname)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output-directory",
                                help="The directory to be used for storing output data. Required.")

    parser.add_argument("-i", "--index-file",
                                help="File containing haloplex indexes (one per line, optional name " \
                                     "in second column separated by tab, comma, or semicolon).")
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
    parser.add_argument("-t", "--top-indexes", type=int,
                                help="Find the n most common indexes. Pair with -l (index length) and -r (read index file). Does not perform any demultiplexing.")
    parser.add_argument("-l", "--index-length", type=int,
                                help="The length of the index.")
    parser.add_argument("-m", "--mismatches", type=int, dest="max_mismatches", default=1,
                                help="The maximum number of mismatches allowed when performing error correction. Default is 1; set to 0 for max speed.")
    parser.add_argument("-p", "--progress-interval", type=int, default=1000,
                                help="Update progress, estimated completion time every N reads (default 1000).")
    arg_vars = vars(parser.parse_args())
    # It's my namespace and I'll clobber it if I want to
    locals().update(arg_vars)
    if arg_vars['top_indexes']:
        if not arg_vars['index_length']:
            raise SyntaxError("Must indicate index length to tally.")
        if not arg_vars['read_index']:
            raise SyntaxError("Must indicate file to parse for indexes.")
        else:
            count_top_indexes(top_indexes, read_index, index_length, progress_interval)
    else:
        main(read_one, read_two, read_index, data_directory, read_index_num, output_directory, index_file, max_mismatches, force_overwrite, progress_interval)
