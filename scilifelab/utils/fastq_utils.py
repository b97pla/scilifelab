"""Utilities for handling FastQ data"""
import gzip
import os
import re
from scilifelab.illumina.hiseq import HiSeqRun
         
class FastQParser:
    """Parser for fastq files, possibly compressed with gzip. 
       Iterates over one record at a time. A record consists 
       of a list with 4 elements corresponding to 1) Header, 
       2) Nucleotide sequence, 3) Optional header, 4) Qualities"""
    
    def __init__(self,file,filter=None):
        self.fname = file
        self.filter = filter
        fh = open(file,"rb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:
            self._fh = fh
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
    """Writes fastq records, where each record is a list with 4 elements
       corresponding to 1) Header, 2) Nucleotide sequence, 3) Optional header, 
       4) Qualities. If the supplied filename ends with .gz, the output file 
       will be compressed with gzip"""
       
    def __init__(self,file):
        self.fname = file
        fh = open(file,"wb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:    
            self._fh = fh
        self._records_written = 0
        
    def name(self):
        return self.fname
    
    def write(self,record):
        self._fh.write("{}\n".format("\n".join([r.strip() for r in record])))
        self._records_written += 1
    
    def rwritten(self):
        return self._records_written
    
    def close(self):
        self._fh.close()

class BarcodeExtractor():
    """Parse a FastQ-file and extract the barcode assumed to be at the 
       given offset and of specified length
    """
    
    def __init__(self,  fqfile, casava18=True, offset=101, length=6):
        fh = open(fqfile)
        if os.path.splitext(fqfile)[1] == ".gz":
            self.fh = gzip.GzipFile(fileobj=fh)
        else:
            self.fh = fh
        self.start = offset
        self.end = offset+length
        self.casava18 = casava18
        self._next = self.setup_next(self.start, self.end)
        
    def __iter__(self):
        return self
    def next(self):
        return self._next(self)
    
    def setup_next(self, start, end):
        """Return the function to extract the barcode
        """
        if not self.casava18:
            def _next(self):
                _, seq, _, _ = [self.fh.next() for i in xrange(4)]
                return seq[start:end]
        else:
            def _next(self):
                header, _, _, _ = [self.fh.next() for i in xrange(4)]
                return header.strip().rsplit(":",1)[1]
        return _next


def avgQ(record,offset=33):
    qual = record[3].strip()
    l = len(qual)
    s = sum([ord(c) for c in qual])
    return round(float(s - l*offset)/l,1)
    
def gtQ30(record,offset=33):
    qual = record[3].strip()
    cutoff = 30 + offset
    g = 0
    for c in qual:
        if ord(c) >= cutoff:
            g += 1
    return round(100*float(g)/len(qual),1)

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

def is_read_pair(rec1, rec2, casava18=True):
    """Returns true if the two records belong to the same read pair, determined by matching the header strings and disregarding
       the read field
    """
    # Handle pre-casava1.8 headers
    if not casava18:
        return (rec1[0][0:-1] == rec2[0][0:-1])
    
    r1 = rec1[0].split(' ')
    r2 = rec2[0].split(' ')
    return (len(r1) == 2 and len(r2) == 2 and r1[0] == r2[0] and r1[1][1:] == r2[1][1:])

def demultiplex_fastq(outdir, samplesheet, fastq1, fastq2=None):
    """Demultiplex a bcl-converted illumina fastq file. Assumes it has the index sequence
    in the header a la CASAVA 1.8+
    """
    outfiles = {}
    counts = {}
    sdata = HiSeqRun.parse_samplesheet(samplesheet)
    reads = [1]
    if fastq2 is not None:
        reads.append(2)
        
    # For each Lane-Index combination, create a file and open a filehandle
    for sd in sdata:
        lane = sd['Lane']
        index = sd['Index']
        if lane not in outfiles:
            outfiles[lane] = {}
            counts[lane] = {}
        outfiles[lane][index] = []
        counts[lane][index] = 0
        for read in reads:
            fname = "tmp_{}_{}_L00{}_R{}_001.fastq.gz".format(sd['SampleID'],
                                                              index,
                                                              lane,
                                                              read)
            outfiles[lane][index].append(FastQWriter(os.path.join(outdir,fname)))
    
    # Parse the input file(s) and write the records to the appropriate output files
    fhs = [FastQParser(fastq1)]
    if fastq2 is not None:
        fhs.append(FastQParser(fastq2))
        
    for r, fh in enumerate(fhs):
        for record in fh:
            header = parse_header(record[0])
            lane = str(header['lane'])
            index = header['index']
            if lane in outfiles and index in outfiles[lane]:
                outfiles[lane][index][r].write(record)
                counts[lane][index] += 1
    
    # Close filehandles and replace the handles with the file names
    for lane in outfiles.keys():
        for index in outfiles[lane].keys():
            for r, fh in enumerate(outfiles[lane][index]):
                fh.close()
                fname = fh.name()
                # If no sequences were written, remove the temporary file and the entry from the results
                if counts[lane][index] == 0:
                    os.unlink(fname)
                    del outfiles[lane][index]
                    break
                
                # Rename the temporary file to a persistent name
                nname = fname.replace("tmp_","")
                os.rename(fname,nname)
                outfiles[lane][index][r] = nname
    
    return outfiles

  
def create_final_name(fname, date, fc_id, sample_name):
    """Create the final name of the delivered file
    """
    
     # Split the file name according to CASAVA convention
    m = re.match(r'(\S+?)_(?:[ACGTN\-]+|NoIndex|Undetermined)_L0*(\d+)_R(\d)_\d+\.fastq(.*)', fname)
    if m is not None:
        lane = m.group(2)
        read = m.group(3)
        ext = m.group(4)
    else:
        # Split the file name according to bcbb convention
        m = re.match(r'(\d+)_(\d+)_([^_]+)_(\d+)_(?:nophix_)?(\d+)_fastq.txt(.*)', fname)
        if m is None:
            raise ValueError("Could not parse file name {:s} correctly!".format(fname))
        lane = m.group(1)
        read = m.group(5)
        ext = m.group(6)
            
    dest_file_name = "{:s}.fastq{:s}".format("_".join([lane,
                                                       date,
                                                       fc_id,
                                                       sample_name,
                                                       read]),
                                             ext.replace('..','.'))
    return dest_file_name


    
    
    
