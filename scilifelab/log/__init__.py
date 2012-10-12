"""
log module
"""
from logbook import Logger

LOG = Logger('scilifelab')

def minimal_logger(name):
    return Logger(name)
