"""
Flowcell object
"""

import os
import sys
import yaml
import json
import csv
import glob
from cStringIO import StringIO

from pmtools.lib.runinfo import _yaml_to_tab
    
class Flowcell(object):
    lanelabels = ["lane", "description", "flowcell_id", "genome_build", "analysis"]
    mplabels = ["barcode_type", "barcode_id", "sample_prj", "name", "sequence"]
    header = lanelabels + mplabels
    runinfo_tab = None
    
    def __init__(self, infile=None):
        if not infile:
            return
        with open(infile) as fh:
            runinfo_yaml = yaml.load(fh)
        self.runinfo_tab = self._yaml_to_tab(runinfo_yaml)

    def __repr__(self):
        return str(self.runinfo_tab)
    
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
        for row in self.runinfo_tab:
            (lane, description, flowcell_id, genome_build, analysis,
             barcode_type, barcode_id, sample_prj, name, sequence) = row
            if not yaml_out.has_key(lane):
                yaml_out[lane] = dict(lane=lane, genome_build=genome_build, description=description,
                                      analysis=analysis, flowcell_id=flowcell_id, multiplex=[])
            yaml_out[lane]["multiplex"].append(dict(barcode_id=barcode_id, barcode_type=barcode_type, genome_build=genome_build,
                                                    name=name, sample_prj=sample_prj, sequence=sequence))
        return yaml.dump(yaml_out.values())
        
    def _get_rows(i):
        return [row for row in self.runinfo_tab[i:]]

    def _column(label):
        i = self.header.index(label)
        return [row[i] for row in self.runinfo_tab]

    def projects(self):
        """List flowcell projects"""
        return list(set(self._column("sample_prj")))

    def lanes(self):
        """List flowcell lanes"""
        return list(set(self._column(runinfo_tab, "lane")))

    def barcodes(self):
        """List barcodes"""
        pass

    def subset(self, column, query):
        """Subset runinfo. Returns new flowcell object."""
        pruned_fc = Flowcell()
        vals = list(self._column(column))
        i = [j + 1 for j in range(0, len(vals)) if vals[j]==query]
        pruned_fc.runinfo_tab = [runinfo_tab[j] for j in i]
        return pruned_fc

