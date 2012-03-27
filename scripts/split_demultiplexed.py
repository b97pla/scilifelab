"""
Script that will split fastq files with demultiplex information in the (CASAVA 1.8+ formatted) header
"""

import os
import csv
import re
import operator
from miseq_data import (FastQParser, FastQWriter, MiSeqSampleSheet, group_fastq_files)
from optparse import OptionParser

def main(fastq_files, outdir, samplesheet):
    
    samples = {}
    if samplesheet:
        ss = MiSeqSampleSheet(samplesheet)
        names = ss.sample_names()
        names.insert(0,"unmatched")
        for i,name in enumerate(names):
            samples[str(i)] = name
            
    _split_fastq_batches(group_fastq_files(fastq_files),outdir,samples)
        
def _split_fastq_batches(inputs, outdir, samples={}):
            
    # Loop over the fastq files
    for fastq_files in inputs:
        fastq_names = [os.path.basename(f) for f in fastq_files]
        prefix = os.path.commonprefix(fastq_names).strip("_")
        suffix = os.path.commonprefix([f[::-1] for f in fastq_names])[::-1]
            
        counts = _split_fastq(fastq_files,outdir,prefix,suffix,samples)
            
    # Write the multiplex metrics
    prefix = os.path.commonprefix([os.path.basename(f) for f in reduce(operator.add,inputs)]).strip("_")
    metrics_file = _write_metrics(counts,outdir,prefix,samples)
    
def _split_fastq(fastq_input, outdir, outprefix, outsuffix, samples):

    if not os.path.exists(outdir):
        os.mkdir(outdir) 
    
    out_handles = {}    
    for file in fastq_input:
        iter = FastQParser(file)
        for record in iter:
            index = record[0].rfind(":")
            i = record[0][index+1:].strip()
            # open a file handle to the index file if it's not already available
            if i not in out_handles:
                out_file = os.path.join(outdir,"%s_%s%s" % (outprefix,samples.get(i,i),outsuffix))
                out_handles[i] = FastQWriter(out_file)
            out_handles[i].write(record)
    
    # summarize the written records and close the file handles
    counts = {}
    for i,oh in out_handles.items():
        counts[i] = oh.rwritten()
        oh.close()
        
    return counts
    
def _write_metrics(counts, outdir, outprefix, samples):
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    metrics_file = os.path.join(outdir,"%s.demultiplex_metrics" % outprefix)
    with open(metrics_file,"wb") as fh:
        cw = csv.writer(fh,dialect=csv.excel_tab)
        cw.writerow(["index","samplesheet sample name","records"])
        for index, count in counts.items():
            cw.writerow([index,samples.get(index,"N/A"),count])

    return metrics_file
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", "--outdir", dest="outdir", default=os.getcwd())
    parser.add_option("-s", "--samplesheet", dest="samplesheet", default={})
    options, args = parser.parse_args()
    
    main(args,options.outdir,options.samplesheet)
