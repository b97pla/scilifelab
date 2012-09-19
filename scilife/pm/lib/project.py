"""
project library

Provides functionality for project folder structure

Project folders should adhere to a common format. Currently
the format is

project
  project_git
    config
      tool-data          [optional]
    sbatch
    doc
  biodata                [optional]
  data
    flowcell_id1
    flowcell_id2
    ...
  intermediate
    analysis_id1
    analysis_id2
    ...

For historical reasons, data and intermediate may be prefixed with nobackup.
"""
