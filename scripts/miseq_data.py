#!/usr/bin/env python
import os
import sys
from optparse import OptionParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip



def main(samplesheet):
    ssobj = MiSeqSampleSheet(samplesheet)

def count_reads(readfiles):
    
    counts = {}
    for readfile in readfiles:
        with gzip.open(readfile) as fh:
            for (name,_,_) in FastqGeneralIterator(fh):
                index = name.rfind(":")
                sample = int(name[index+1:].strip())
                if sample not in counts:
                    counts[sample] = 0
                counts[sample] += 1
    total = 0
    for sample, count in counts.items():
        total += count
    counts["total"] = total
    return counts

class MiSeqSampleSheet:
    
    def __init__(self, ss_file):
        assert os.path.exists(ss_file), "Samplesheet %s does not exist" % ss_file
        
        print count_reads([ss_file])
        sys.exit(0)
        setattr(self,"samplesheet",ss_file)
        self._parse_sample_sheet()
        
    def _parse_sample_sheet(self):
        
        # Parse the samplesheet file into a data structure
        data = {}
        with open(self.samplesheet,"r") as fh:
            current = None
            for line in fh:
                line = line.strip()
                if line.startswith("["):
                    current = line.strip("[], ")
                    data[current] = {}
                else:
                    [opt,val] = line.split(",",1)
                    data[current][opt] = val
    
        # Assign the parsed attributes to class attributes
        for option, value in data.get("Header",{}).items():
            setattr(self, option, value)
        for option, value in data.get("Settings",{}).items():
            setattr(self, option, value)
        
        # Parse sample data
        first_data_col = "Sample_ID"
        if "Data" in data and first_data_col in data["Data"]:
            data_header = data["Data"][first_data_col].split(",")
            samples = {}
            for sample_id, sample_data in data["Data"].items():
                if sample_id == first_data_col: continue
                samples[sample_id] = dict(zip(data_header,sample_data.split(",")))
                samples[sample_id][first_data_col] = sample_id
            setattr(self, "samples", samples)

    def sample_field(self, sample_id, sample_field=None):
        samples = getattr(self,"samples",{})
        assert sample_id in samples, "The sample '%s' was not found in samplesheet %s" % (sample_id,self.samplesheet)
        if sample_field is None:
            return samples[sample_id]
        assert sample_field in samples[sample_id], "The sample field '%s' was not found in samplesheet %s" % (sample_field,self.samplesheet)
        return samples[sample_id][sample_field] 

    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--samplesheet", dest="samplesheet", default=None)
    #parser.add_option("-n", "--dry-run", dest="dryrun", action="store_true", default=False)
    options, args = parser.parse_args()
    
    
    main(options.samplesheet)
