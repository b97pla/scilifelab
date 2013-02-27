""" A module for handling MiSeq-specific files and folders
"""
import os
import re
import glob
from collections import OrderedDict
from scilifelab.bcbio.qc import FlowcellRunMetricsParser
from scilifelab.illumina.hiseq import HiSeqSampleSheet
    
def group_fastq_files(fastq_files):
    """Divide the input fastq files into batches based on lane and read, ignoring set"""
        
    regexp = r'_(L\d+)_([RI]\d+)_'
    batches = {}
    for fastq_file in fastq_files:
        m = re.search(regexp, fastq_file)
        if not m or len(m.groups()) < 2:
            print "WARNING: Could not determine lane and read from input file %s" % fastq_file
            continue
        
        batch = "%s%s" % (m.group(1).strip(),m.group(2).strip())
        if batch not in batches:
            batches[batch] = []
        batches[batch].append(fastq_file)

    return batches.values()


class MiSeqRun:
    def __init__(self, run_dir):
        self._run_dir = os.path.normpath(run_dir)
        assert os.path.exists(self._run_dir), "The path %s is invalid" % self._run_dir
        ss_file = self._find_samplesheet()
        if ss_file is not None:
            samplesheet = MiSeqSampleSheet(ss_file)
            self.samplesheet = samplesheet
        
        parser = FlowcellRunMetricsParser(self._run_dir)
        self.run_config = parser.parseRunParameters()
        self._fastq = self._fastq_files()
        
    def write_hiseq_samplesheet(self, samplesheet):
        """Export the metadata for this run in a HiSeq samplesheet format
        """
        hs_ssheet = HiSeqSampleSheet(self.samplesheet.to_hiseq(self.run_config))
        hs_ssheet.write(samplesheet)
        
    def _data_dir(self):
        return os.path.join(self._run_dir,"Data")
    def _intensities_dir(self):
        return os.path.join(self._data_dir(),"Intensities")
    def _basecalls_dir(self):
        return os.path.join(self._intensities_dir(),"BaseCalls")
    def _multiplex_dir(self):
        return os.path.join(self._basecalls_dir(),"Multiplex")
    def _alignment_dir(self):
        return os.path.join(self._basecalls_dir(),"Alignment")
    def _runParameters(self):
        return os.path.join(self._run_dir,"runParameters.xml")
    
    def _fastq_files(self, fastq_dir=None):
        if fastq_dir is None:
            fastq_dir = self._basecalls_dir()
        
        fastq_files = group_fastq_files(glob.glob(os.path.join(fastq_dir,"*.fastq*")))
        return fastq_files
    
    def _find_samplesheet(self):
        dirs = [self._run_dir,
                self._basecalls_dir()]
        for dir in dirs:
            ss = os.path.join(dir,"SampleSheet.csv")
            if os.path.exists(ss):
                return ss
        return None
    
    def _split_fastq(self):
        
        samples = self.samplesheet.sample_names()
        samples.insert(0,"unmatched")
        sample_names = {}
        for i,name in enumerate(samples):
            sample_names[str(i)] = name
        
        out_dir = self._multiplex_dir()
        
        import split_demultiplexed 
        split_demultiplexed._split_fastq_batches(self._fastq,out_dir,sample_names)


class MiSeqSampleSheet:
    def __init__(self, ss_file):
        assert os.path.exists(ss_file), \
            "Samplesheet %s does not exist" % ss_file

        setattr(self, "samplesheet", ss_file)
        self.data_header = ["Sample_ID",
                            "Sample_Name",
                            "Sample_Plate",
                            "Sample_Well",
                            "Sample_Project",
                            "index",
                            "I7_Index_ID",
                            "index2",
                            "I5_Index_ID",
                            "Description",
                            "Manifest",
                            "GenomeFolder"]
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
                    if current is None:
                        current = "NoSection"
                    s = line.split(",",1)
                    if len(s) > 1: 
                        data[current][s[0]] = s[1]
                    else:
                        data[current][line] = ''
    
        # Assign the parsed attributes to class attributes
        for option, value in data.get("Header",{}).items():
            setattr(self, option.replace(" ", ""), value)

        for option, value in data.get("Settings",{}).items():
            setattr(self, option, value)
        if "Data" not in data:
            data["Data"] = {}
            data["Data"][self.data_header[0]] = ",".join(self.data_header[1:])
            for option, value in data.get("NoSection",{}).items():
                data["Data"][option] = value
            
        # Parse sample data
        first_data_col = "Sample_ID"
        if "Data" in data and first_data_col in data["Data"]:
            self.data_header = [s.lower() for s in data["Data"][first_data_col].split(",")]
            samples = {}
            for sample_id, sample_data in data["Data"].items():
                if sample_id == first_data_col:
                    continue

                samples[sample_id] = dict(zip(self.data_header,sample_data.split(",")))
                samples[sample_id][first_data_col.lower()] = sample_id
                
            setattr(self, "samples", samples)

    def sample_names(self):
        """Return the name of the samples in the same order as they are listed in
        the samplesheet.
        """
        samples = getattr(self,"samples",{})
        
        if getattr(self, "_sample_names", None) is None:
            sample_names = []
            with open(self.samplesheet,"r") as fh:
                for line in fh:
                    if line.startswith("[Data]"):
                        for line in fh:
                            data = line.split(",")
                            if len(data) == 0 or data[0].startswith("["):
                                break

                            if data[0] in samples:
                                sample_names.append(data[0])

            self._sample_names = sample_names
        
        return self._sample_names
        
        
    def sample_field(self, sample_id, sample_field=None):
        samples = getattr(self,"samples",{})
        assert sample_id in samples, \
            "The sample '%s' was not found in samplesheet %s" % (sample_id,self.samplesheet)
        if sample_field is None:
            return samples[sample_id]

        assert sample_field in samples[sample_id], \
            "The sample field '%s' was not found in samplesheet %s" % (sample_field,self.samplesheet)

        return samples[sample_id][sample_field]

    def to_hiseq(self, run_config={}):
        """Convert Miseq SampleSheet to HiSeq formatted Samplesheet.
        """
        FCID = run_config.get('Barcode','NA')
        Lane = "1"
        SampleRef = "NA"
        Description = "NA"
        Control = "N"
        Recipe = "NA"
        Operator = "NA"

        rows = []
        for sampleID, info in self.samples.iteritems():
            row = OrderedDict()
            row["FCID"] = FCID
            row["Lane"] = Lane
            row["SampleID"] = sampleID
            row["SampleRef"] = self._extract_reference_from_path(info.get('genomefolder',''))
            row["Index"] = info.get('index','NoIndex')
            row["Description"] = info.get('description','')
            row["Control"] = Control
            row["Recipe"] = Recipe
            row["Operator"] = Operator
            row["SampleProject"] = info.get('sample_project','Unknown')

            rows.append(row)

        return rows
    
    def _extract_reference_from_path(self, path):
        """Attempts to extract a name of a reference assembly from a path
        """
    
        head = path
        regexp = r'[a-zA-Z]+[0-9\.]+$'
        while head is not None and len(head) > 0:
            head, tail = os.path.split(head.replace('\\','/'))
            if re.match(regexp, tail) is not None:
                return tail
            
        return path
    
    