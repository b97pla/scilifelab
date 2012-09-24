"""
Flowcell object
"""

import os
import sys
import re
import yaml
import json
import csv
import glob
import copy
from cStringIO import StringIO
from scilifelab.pm.utils.misc import filtered_walk

## FIX ME: what should be returned from object functions, and what
## should be done behind the scenes?

## FIX ME: make generic flowcell object, that then Illumina, MiSeq,
## SOLiD subclass from

class Flowcell(object):
    """Class for handling (Illumina) run information.

    Run information is stored as a table."""
    ## Run information
    fc_name = None
    fc_date = None
    unique_lanes = False
    ## Unique column names
    _keys = dict(lane = ['lane', 'lane_description', 'flowcell_id', 'lane_analysis', 'genome_build'],
                 mp = ['mp_analysis', 'barcode_id', 'barcode_type', 'sample_prj', 'name', 'sequence', 'files', 'genomes_filter_out', 'mp_description', 'results'])
    keys = _keys['lane'] + _keys['mp']
    
    ## csv keys
    # _csv_keys = dict('FCID'='flowcell_id', 'Lane' = 'lane', 'SampleID'='name', 'SampleRef' = 'genome_build', 'Index'='sequence', 'Description' = 'sample_prj',
    #                  'Control'=None, 'Recipe'=None, 'Operator'=None, 'SampleProject'='sample_prj')
    # _csv_keys = ['FCID', 'Lane', 'SampleID', 'SampleRef', 'Index', 'Description', 'Control', 'Recipe', 'Operator', 'SampleProject']
    _csv_keys = ['flowcell_id', 'lane', 'name', 'genome_build', 'sequence', 'sample_prj', 'control', 'recipe', 'operator', 'sample_prj']

    ## sample keys
    samples = dict()
    ## lane files
    lane_files = dict()
    
    # keys to be printed for yaml output
    _out_yaml_keys = dict(lane= ['lane', 'lane_description', 'flowcell_id', 'lane_analysis', 'genome_build'],
                     mp = ['mp_analysis', 'barcode_id', 'barcode_type', 'sample_prj', 'name', 'sequence', 'files', 'genomes_filter_out', 'mp_description'])
    ## keys to write for table
    _out_table_columns = ['lane', 'lane_description', 'flowcell_id', 'lane_analysis', 'genome_build', 'barcode_id', 'barcode_type', 'sample_prj', 'name', 'sequence', 'genomes_filter_out']

    def __init__(self, infile=None):
        self.filename = None
        self.path = None
        self.data = None
        self.i = 0
        if not infile:
            return
        self.data = self._read(infile)
        self._set_sample_dict()

    def fc_id(self):
        m = re.search("([0-9]+)_[A-Za-z0-9]+_[A-Za-z0-9]+_([A-Z0-9]+)", os.path.dirname(self.filename))
        if m:
            (self.fc_date, self.fc_name) = (m.group(1), m.group(2))
        return "{}_{}".format(self.fc_date, self.fc_name)
        
    def __iter__(self):
        return self

    def next(self):
        if self.i >= len(self.data):
            self.i = 0
            raise StopIteration
        row = self.data[self.i]
        self.i = self.i + 1 
        return dict(zip(self.keys, row))

    def __repr__(self):
        return "Flowcell(filename={})".format(self.filename)

    def __str__(self):
        fh = StringIO()
        w = csv.writer(fh, delimiter="\t", quoting=True)
        tab_out = []
        for row in self.data:
            h = self._keys['lane'] + self._keys['mp']
            d = dict(zip(h, row))
            if not d.get("barcode_id", None):
                continue
            tab_out.append([d[k] for k in self._out_table_columns])
        w.writerows([self._out_table_columns] + tab_out)
        return fh.getvalue()

    def __len__(self):
        if not self.data:
            return 0
        return len(self.data)

    def _read(self, infile):
        """Read infile. Pass to correct read wrapper."""
        (root, ext) = os.path.splitext(infile)
        if ext in [".csv"]:
            return self._read_csv(infile)
        else:
            return self._read_yaml(infile)

    def _read_csv(self, infile):
        """Reads samplesheet csv file. Returns self converted to tab"""
        out = []
        if not os.path.exists(infile):
            return None
        first = True
        with open(infile) as fh:
            runinfo_csv = csv.reader(fh)
            for row in runinfo_csv:
                if first:
                    first = False
                    continue
                d = dict(zip(self._csv_keys, row))
                d['analysis'] = "Align_standard_seqcap"
                d['lane_description'] = "Lane {}, {}".format(d['lane'], d['sample_prj'].replace("__", "."))
                d['lane_analysis'] = None
                d['mp_description'] = "{}_{}".format(d['sample_prj'].replace("__", "."), d['name'])
                d['mp_analysis'] = None
                d['barcode_id'] = None
                d['barcode_type'] = "Samplesheet"
                d['files'] = [] ##["{}_{}_L00{}_R1_001.fastq".format(d['name'], d['sequence'], d['lane']), "{}_{}_L00{}_R2_001.fastq".format(d['name'], d['sequence'], d['lane'])]
                d['genomes_filter_out'] = None ##'phix'
                d['results'] = None
                newrow = [d[k] for k in self.keys]
                out.append(newrow)
        self.filename = os.path.abspath(infile)
        return out
        
    def _read_yaml(self, infile):
        """Reads runinfo yaml file. Returns self converted to tab"""
        if not os.path.exists(infile):
            return None
        with open(infile) as fh:
            runinfo_yaml = yaml.load(fh)
        self.filename = os.path.abspath(infile)
        return self._yaml_to_tab(runinfo_yaml)
    
    def as_yaml(self):
        return self._tab_to_yaml()

    def _set_sample_dict(self):
        i = 0
        self.samples = {}
        for row in self.data:
            d = dict(zip(self.keys, row))
            key = "{}_{}".format(d['lane'], d['barcode_id'])
            self.samples[key] = i
            i += 1
            
    def _yaml_to_tab(self, runinfo_yaml):
        """Convert yaml to internal representation"""
        out = []
        for info in runinfo_yaml:
            self.lane_files[info.get('lane', None)] = []
            laneinfo = [info.get(x.replace("lane_", "")) for x in self._keys['lane']]
            for mp in info.get("multiplex", None):
                mpinfo = [mp.get(x.replace("mp_", "")) for x in self._keys['mp']]
                line = laneinfo + mpinfo
                line[self.keys.index("files")] = []
                line[self.keys.index("results")] = []
                out.append(line)
        return out
    
    def _tab_to_yaml(self):
        """Convert internal representation to yaml"""
        yaml_out = dict()
        for row in self.data:
            d = dict(zip(self.keys, row))
            if not d.get("barcode_id", None):
                continue
            if not yaml_out.has_key(d['lane']):
                yaml_out[d['lane']] = dict((k.replace("lane_", ""), d[k]) for k in self._out_yaml_keys['lane'] if k in d and not d[k] is None)
                if not self.unique_lanes:
                    yaml_out[d['lane']]["multiplex"] = []
            d_mp = dict((k.replace("mp_", ""), d[k]) for k in self._out_yaml_keys['mp'] if k in d and not d[k] is None)
            ## Fix description, results, files and analysis
            d_mp["analysis"] = d_mp.get("analysis", yaml_out[d['lane']].get("analysis", None))
            d_mp["description"] = d_mp.get("description", "{}_{}".format(str(d_mp.get("sample_prj", None)), str(d_mp.get("name", None))))
            if d_mp.get("files", None):
                d_mp["files"] = list(set(d_mp["files"]))
            if not self.unique_lanes:
                yaml_out[d['lane']]["multiplex"].append(d_mp)
            else:
                yaml_out[d['lane']]['description'] = d_mp['description']
        yaml_out_final = dict(details=yaml_out.values())
        if self.fc_name:
            yaml_out_final['fc_name'] = self.fc_name
        if self.fc_date:
            yaml_out_final['fc_date'] = self.fc_date
        return yaml.dump(yaml_out_final)
    
    def get_sample(self, key):
        return self.data[self.samples[key]]

    def get_entry(self, key, label):
        return self.data[self.samples[key]][self.keys.index(label)]

    def set_entry(self, key, label, value):
        self.data[self.samples[key]][self.keys.index(label)] = value

    def append_to_entry(self, key, label, value):
        if not self.data[self.samples[key]][self.keys.index(label)]:
            self.data[self.samples[key]][self.keys.index(label)] = []
        self.data[self.samples[key]][self.keys.index(label)].append(value)
                
    def _column(self, label):
        i = self.keys.index(label)
        return [row[i] for row in self.data if not row[i] is None]

    def projects(self):
        """List flowcell projects"""
        return list(set(self._column("sample_prj")))

    def lanes(self):
        """List flowcell lanes"""
        return list(set(self._column("lane")))

    def barcodes(self, lane):
        """List barcodes for a lane"""
        fc = self.subset("lane", lane)
        return list(fc._column("barcode_id"))

    def names(self, lane):
        """List names for a lane"""
        fc = self.subset("lane", lane)
        return list(fc._column("name"))

    def barcode_id_to_name(self, lane):
        """Map barcode id to name"""
        return dict(zip(self.barcodes(lane), self.names(lane)))

    def barcode_name_to_id(self, lane):
        """Map barcode name to id"""
        return dict(zip(self.names(lane), self.barcodes(lane)))

    def fc_with_unique_lanes(self):
        """Transform flowcell to one with unique lane numbers"""
        new_fc = copy.deepcopy(self)
        new_fc.filename = self.filename.replace(".yaml", "-unique-lane.yaml")
        lane_index = new_fc.keys.index("lane")
        lane = 1
        for j in range(0, len(new_fc.data)):
            new_fc.data[j][lane_index] = str(lane)
            lane = lane + 1
        new_fc.unique_lanes = True
        return new_fc
            
    def subset(self, column, query):
        """Subset runinfo. Returns new flowcell object."""
        pruned_fc = Flowcell()
        vals = list(self._column(column))
        i = [j for j in range(0, len(vals)) if vals[j]==query]
        pruned_fc.data = [self.data[j] for j in i]
        pruned_fc.filename = self.filename.replace(".yaml", "-pruned.yaml")
        pruned_fc._set_sample_dict()
        pruned_fc.lane_files = dict((x, self.lane_files[x]) for x in pruned_fc.lanes())
        return pruned_fc

    def load(self, paths, runinfo="run_info.yaml"):
        """Load run information.

        :param: paths - one or several paths to look in
        """
        for p in paths:
            if not os.path.exists(p):
                next
            data = self._read(os.path.join(p, runinfo))
            if not data is None:
                self.data = data
                self.filename = os.path.join(p, runinfo)
                self._set_sample_dict()
                break
        if not data:
            return None
        else:
            return data

    def glob_pfx_dict(self, ext="", sample=True):
        """Return glob sample prefix regular expression strings as a dict"""
        glob_pfx = dict()
        for smp in self:
            if sample:
                re_str = "{}_{}".format(smp['lane'], smp['barcode_id'])
            else:
                re_str = "{}".format(smp['lane'])
            pattern = "{}_[0-9]+_.?{}(_nophix)?_{}*{}".format(smp['lane'], smp['flowcell_id'], smp['barcode_id'], ext)
            glob_pfx[sample] = pattern
        return glob_pfx

    def glob_pfx_str(self, ext=""):
        """Return glob prefix regular expression strings"""
        glob_pfx = []
        for sample in self:
            pattern = "{}_[0-9]+_.?{}(_nophix)?_{}*{}".format(sample['lane'], sample['flowcell_id'], sample['barcode_id'], ext)
            glob_pfx.append(pattern)
        return glob_pfx
        
    def glob_pfx_re(self, ext=""):
        """Return glob prefix regular expressions."""
        glob_pfx = []
        for sample in self:
            pattern = re.compile("{}_[0-9]+_.?{}(_nophix)?_{}-*{}".format(sample['lane'], sample['flowcell_id'], sample['barcode_id'], ext))
            glob_pfx.append(pattern)
        return glob_pfx

    ## Lane objects
    # 1_120829_AA001AAAXX_nophix.bc_metrics
    # 1_120829_AA001AAAXX_nophix.filter_metrics
    # 1_120829_AA001AAAXX_nophix_1_fastq.txt
    ## Sample files
    # 1_120829_AA001AAAXX_nophix_10.bam
    # 1_120829_AA001AAAXX_nophix_10.sam
    # 1_120829_AA001AAAXX_nophix_10_1.sai
    # 1_120829_AA001AAAXX_nophix_10_2.sai
    # 1_120829_AA001AAAXX_nophix_10_1_fastq.txt
    # 1_120829_AA001AAAXX_nophix_1-sort-dup.align_metrics
    # 1_120829_AA001AAAXX_nophix_1-sort-dup.bam	
    # 1_120829_AA001AAAXX_nophix_1-sort-dup.dup_metrics
    # 1_120829_AA001AAAXX_nophix_1-sort-dup.hs_metrics
    # 1_120829_AA001AAAXX_nophix_1-sort-dup.insert_metrics
    # 1_120829_AA001AAAXX_nophix_1-sort.bam
    def classify_file(self, f):
        re_lane = re.compile('^([0-9]+)_[0-9]+_[A-Za-z0-9]+(_nophix)?\.(filter|bc)_metrics|^([0-9]+)_[0-9]+_[A-Za-z0-9]+(_nophix)?_[12]_fastq.txt')
        m_lane = re_lane.search(os.path.basename(f))
        if m_lane:
            lane = os.path.basename(f).split("_")[0]
            sample = None
            self.lane_files[lane].append(os.path.abspath(f))
            return
        re_sample = re.compile('^([0-9]+)_[0-9]+_[A-Za-z0-9]+(_nophix)?_([0-9]+|unmatched).*')
        m_sample = re_sample.search(os.path.basename(f))
        if m_sample:
            lane = m_sample.group(1)
            sample = m_sample.group(3)
            key = "{}_{}".format(lane, sample)
            if sample == "unmatched":
                self.lane_files[lane].append(os.path.abspath(f))
                return
            if f.find("fastq.txt") > 0:
                self.append_to_entry(key, "files", os.path.abspath(f))
            else:
                self.append_to_entry(key, "results", os.path.abspath(f))
        else:
            return 
        return

    def collect_files(self, path, project=None):
        """Collect files for a given project"""
        if project:
            fc = self.subset("sample_prj", project)
        else:
            fc = self
        pattern = "|".join(fc.glob_pfx_str())
        def file_filter(f):
            if not pattern:
                return
            return re.search(pattern, f) != None
        flist = filtered_walk(path, file_filter)
        for f in flist:
            self.classify_file(f)
        fc.path = path
        return fc
