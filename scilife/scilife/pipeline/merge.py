"""Handle multiple samples present on a single flowcell

From bcbio.pipeline.merge:

This provides functionality for cases where you have the same sample in multiple lanes
on a flowcell. They can be combined from barcoded subsets and are identified
based on unique sample names.
"""
import os
import collections

from bcbio import utils

def organize_samples(dirs, fc_name, fc_date, run_items):
    pass
