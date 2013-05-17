"""pm bcbio utils module"""
import re

FLOWCELL_DIRECTORY_PATTERN = "[0-9]+_[0-9A-Za-z]+_[0-9]+_[A-Z]?[A-Z0-9\-]+"

def validate_fc_directory_format(flowcell):
    """Validate that path conforms to bcbio flowcell directory name"""
    m = re.search(FLOWCELL_DIRECTORY_PATTERN, flowcell)
    if not m:
        return False
    else:
        return True

def fc_id(flowcell):
    """Return fc id"""
    pattern = "[0-9]+_[0-9A-Za-z]+_[0-9]+_[A-Z]?([A-Z0-9\-]+)"
    m = re.search(pattern, flowcell)
    if not m:
        return None
    else:
        return m.group(1)

def fc_fullname(flowcell):
    """Return fc name (fc_date_fc_name)"""
    pattern = "([0-9]+)_[0-9A-Za-z]+_[0-9]+_([A-Z0-9\-]+)"
    m = re.search(pattern, flowcell)
    if not m:
        return None
    else:
        return "{}_{}".format(m.group(1), m.group(2))
    
def fc_parts(flowcell):
    """return fc_name and fc_date"""
    pattern = "([0-9]+)_[0-9A-Za-z]+_[0-9]+_([A-Z0-9\-]+)"
    m = re.search(pattern, flowcell)
    if not m:
        return None
    else:
        return (m.group(1), m.group(2))



