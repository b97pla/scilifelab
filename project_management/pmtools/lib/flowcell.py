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


## FIX ME: what should be returned from object functions, and what
## should be done behind the scenes?

## FIX ME: make generic flowcell object, that then Illumina, MiSeq,
## SOLiD subclass from

class Flowcell(object):
    """Class for handling (Illumina) flowcell info"""
    lanelabels = ["lane", "description", "flowcell_id", "genome_build", "analysis"]
    mplabels = ["barcode_type", "barcode_id", "sample_prj", "name", "sequence"]
    header = lanelabels + mplabels
    data = None
    
    def __init__(self, infile=None):
        if not infile:
            return None
        self.data = self._read(infile)

    def __str__(self):
        return str(self.data)

    def __len__(self):
        return len(self.data)

    def _read(self, infile):
        """Reads runinfo. Returns self"""
        if not os.path.exists(infile):
            return None
        with open(infile) as fh:
            runinfo_yaml = yaml.load(fh)
        return self._yaml_to_tab(runinfo_yaml)
    
    def as_yaml(self):
        return self._tab_to_yaml()

    def _yaml_to_tab(self, runinfo_yaml):
        """Convert yaml to internal representation"""
        out = []
        for info in runinfo_yaml:
            laneinfo = [info.get(x, None) for x in self.lanelabels]
            for mp in info.get("multiplex", None):
                mpinfo = [mp.get(x, None) for x in self.mplabels]
                line = laneinfo + mpinfo
                out.append(line)
        return out

    def _tab_to_yaml(self):
        """Convert internal representation to yaml"""
        yaml_out = dict()
        for row in self.data:
            (lane, description, flowcell_id, genome_build, analysis,
             barcode_type, barcode_id, sample_prj, name, sequence) = row
            if not yaml_out.has_key(lane):
                yaml_out[lane] = dict(lane=lane, genome_build=genome_build, description=description,
                                      analysis=analysis, flowcell_id=flowcell_id, multiplex=[])
            yaml_out[lane]["multiplex"].append(dict(barcode_id=barcode_id, barcode_type=barcode_type, genome_build=genome_build,
                                                    name=name, sample_prj=sample_prj, sequence=sequence))
        return yaml.dump(yaml_out.values())
        
    def _get_rows(self, i):
        return [row for row in self.data[i:]]

    def _column(self, label):
        i = self.header.index(label)
        return [row[i] for row in self.data]

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
        return pruned_fc

    def load(self, paths, runinfo="run_info.yaml"):
        """Load run information.

        :param: paths - one or several paths to look in
        """
        for p in paths:
            data = self._read(os.path.join(p, runinfo))
            if not data is None:
                self.data = data
                break
        if not data:
            return None
        else:
            return data

    def glob_pfx_str(self, ext=""):
        """Return glob prefix regular expression strings"""
        glob_pfx = []
        for row in self.data:
            pattern = "{}_[0-9]+_?{}(_nophix)?_{}-*{}".format(row[0], row[2], row[6], ext)
            glob_pfx.append(pattern)
        return glob_pfx
        
    def glob_pfx_re(self, ext=""):
        """Return glob prefix regular expressions."""
        glob_pfx = []
        for row in self.data:
            pattern = re.compile("{}_[0-9]+_?{}(_nophix)?_{}-*{}".format(row[0], row[2], row[6], ext))
            glob_pfx.append(pattern)
        return glob_pfx

        
