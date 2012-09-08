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
from cStringIO import StringIO
from pmtools.utils.misc import filtered_walk, group_bcbb_files

## FIX ME: what should be returned from object functions, and what
## should be done behind the scenes?

## FIX ME: make generic flowcell object, that then Illumina, MiSeq,
## SOLiD subclass from

# details:
# - analysis: Align_standard
#   description: Lane 6, M.Uhlen_12_05
#   flowcell_id: C118PACXX
#   genome_build: rn4
#   lane: '6'
#   multiplex:
#   - analysis: Align_standard
#     barcode_id: 1
#     barcode_type: SampleSheet
#     description: M.Uhlen_12_05_P281_125F_index2
#     files:
#     - P281_125F_index2_CGATGT_L006_R1_001.fastq
#     - P281_125F_index2_CGATGT_L006_R2_001.fastq
#     genome_build: rn4
#     genomes_filter_out: phix
#     name: P281_125F_index2
#     sample_prj: M.Uhlen_12_05
#     sequence: CGATGT
# fc_date: '120821'
# fc_name: BC118PACXX

class Flowcell(object):
    """Class for handling (Illumina) run information.

    Run information is stored as a table."""
    ## Run information
    fc_name = None
    fc_date = None
    _keys = dict(lane = ['lane', 'lane_description', 'flowcell_id', 'lane_analysis', 'genome_build', 'lane_results'],
                 mp = ['mp_analysis', 'barcode_id', 'barcode_type', 'sample_prj', 'name', 'sequence', 'files', 'genomes_filter_out', 'mp_description', 'mp_results'])
    _labels = dict(lane = ['lane', 'description', 'flowcell_id', 'analysis', 'genome_build', 'results'],
                   mp = ['analysis', 'barcode_id', 'barcode_type', 'sample_prj', 'name', 'sequence', 'files', 'genomes_filter_out', 'description', 'results'])
    header = _labels['lane'] + _labels['mp']

    def __init__(self, infile=None):
        self.filename = None
        self.data = None
        self.data_dict = None
        if not infile:
            return None
        self.data = self._read(infile)

    def __str__(self):
        fh = StringIO()
        w = csv.writer(fh, delimiter="\t", quoting=True)
        w.writerows([self.header] + self.data)
        return fh.getvalue()

    def _init_dict(self):
        d = {"{}_{}".format(row[0], row[6]):{"name":row[8], "files":[]} for row in self.data}#, "{}".format(row[0])}
        d2 = {"{}".format(l):{"name":str(l),"files":[]} for l in self.lanes()}
        d.update(d2)
        self.data_dict = d
        
    def __len__(self):
        if not self.data:
            return 0
        return len(self.data)

    def _read(self, infile):
        """Reads runinfo. Returns self"""
        if not os.path.exists(infile):
            return None
        with open(infile) as fh:
            runinfo_yaml = yaml.load(fh)
        self.filename = os.path.abspath(infile)
        return self._yaml_to_tab(runinfo_yaml)
    
    def as_yaml(self):
        return self._tab_to_yaml()

    def _yaml_to_tab(self, runinfo_yaml):
        """Convert yaml to internal representation"""
        out = []
        for info in runinfo_yaml:
            laneinfo = [info.get(x) for x in self._labels['lane']]
            for mp in info.get("multiplex", None):
                mpinfo = [mp.get(x) for x in self._labels['mp']]
                line = laneinfo + mpinfo
                out.append(line)
        return out
    
    def _tab_to_yaml(self):
        """Convert internal representation to yaml"""
        yaml_out = dict()
        for row in self.data:
            h= self._keys['lane'] + self._keys['mp']
            d= dict(zip(h, row))
            if not yaml_out.has_key(d['lane']):
                yaml_out[d['lane']] = dict((k.replace("lane_", ""), d[k]) for k in self._keys['lane'] if k in d and not d[k] is None)
                yaml_out[d['lane']]["multiplex"] = []
            d_mp = dict((k.replace("mp_", ""), d[k]) for k in self._keys['mp'] if k in d and not d[k] is None)
            ## Fix description and analysis
            d_mp["analysis"] = d_mp.get("analysis", yaml_out[d['lane']]["analysis"])
            d_mp["description"] = d_mp.get("description", "{}_{}".format(str(d_mp.get("sample_prj", None)), str(d_mp.get("name", None))))
            yaml_out[d['lane']]["multiplex"].append(d_mp)
        return yaml.dump(yaml_out.values())
    
    def _get_rows(self, i):
        return [row for row in self.data[i:]]

    def _column(self, label):
        i = self.header.index(label)
        return [row[i] for row in self.data]

    # def set_path(self, path):
    #     """Set flowcell path"""
        

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

    def subset(self, column, query):
        """Subset runinfo. Returns new flowcell object."""
        pruned_fc = Flowcell()
        vals = list(self._column(column))
        i = [j for j in range(0, len(vals)) if vals[j]==query]
        pruned_fc.data = [self.data[j] for j in i]
        #pruned_fc.path = self.path
        pruned_fc.filename = self.filename.replace(".yaml", "-pruned.yaml")
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
                break
        if not data:
            return None
        else:
            return data

    def glob_pfx_dict(self, ext="", sample=True):
        """Return glob sample prefix regular expression strings as a dict"""
        glob_pfx = dict()
        for row in self.data:
            if sample:
                re_str = "{}_{}".format(row[0], row[2])
            else:
                re_str = "{}".format(row[0])
            pattern = "{}_[0-9]+_.?{}(_nophix)?_{}*{}".format(row[0], row[2], row[6], ext)
            glob_pfx[sample] = pattern
        return glob_pfx

    def glob_pfx_str(self, ext=""):
        """Return glob prefix regular expression strings"""
        glob_pfx = []
        for row in self.data:
            pattern = "{}_[0-9]+_.?{}(_nophix)?_{}*{}".format(row[0], row[2], row[6], ext)
            glob_pfx.append(pattern)
        return glob_pfx
        
    def glob_pfx_re(self, ext=""):
        """Return glob prefix regular expressions."""
        glob_pfx = []
        for row in self.data:
            pattern = re.compile("{}_[0-9]+_.?{}(_nophix)?_{}-*{}".format(row[0], row[2], row[6], ext))
            glob_pfx.append(pattern)
        return glob_pfx

    ### FIXME: patterns should be set in bcbb_extension
    def get_files(self, path, ftype="", ext=".bam", project=None, lane=None):
        """Get files with a given pattern"""
        if project:
            fc = self.subset("sample_prj", project)
        elif lane:
            pass
        else:
            fc = self
        pattern = "|".join(fc.glob_pfx_str())
        def file_filter(f):
            if not pattern:
                return
            return re.search(pattern, f) != None
        flist = filtered_walk(path, file_filter)
        return flist

    def collect_files(self, path, project=None):
        """Collect files for a given project"""
        if project:
            fc = self.subset("sample_prj", project)
        else:
            fc = self
        if fc.data_dict is None:
            fc._init_dict()
        pattern = "|".join(fc.glob_pfx_str())
        def file_filter(f):
            if not pattern:
                return
            return re.search(pattern, f) != None
        flist = filtered_walk(path, file_filter)
        [fc.data_dict[group_bcbb_files(x)]["files"].append(x.replace(path, "").lstrip("/")) for x in flist]
        return fc
    

        
