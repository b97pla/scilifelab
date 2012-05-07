"""Utilities for handling FastQ data"""
import gzip

PHRED_OFFSET = 33
         
class FastQParser:
    """Parser for fastq files, possibly compressed with gzip. 
       Iterates over one record at a time. A record consists 
       of a list with 4 elements corresponding to 1) Header, 
       2) Nucleotide sequence, 3) Optional header, 4) Qualities"""
    
    def __init__(self,file):
        fh = open(file,"rb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:
            self._fh = fh
        self._records_read = 0
        
    def __iter__(self):
        return self
    def next(self):
        record = []
        for i in range(4):
            record.append(self._fh.next().strip())
        self._records_read += 1
        
        return record

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
        fh = open(file,"wb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:    
            self._fh = fh
        self._records_written = 0
        
    def write(self,record):
        for row in record:
            self._fh.write("%s\n" % row.strip("\n"))
        self._records_written += 1
    
    def rwritten(self):
        return self._records_written
    
    def close(self):
        self._fh.close()

def avgQ(record,offset=PHRED_OFFSET):
    qual = record[3].strip()
    l = len(qual)
    s = sum([ord(c) for c in qual])
    return round(float(s - l*offset)/l,1)
    
def gtQ30(record,offset=PHRED_OFFSET):
    qual = record[3].strip()
    cutoff = 30 + offset
    g = 0
    for c in qual:
        if ord(c) >= cutoff:
            g += 1
    return round(100*float(g)/len(qual),1)
