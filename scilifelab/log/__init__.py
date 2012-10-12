"""
log module
"""
from logbook import Logger

log = Logger('scilifelab')

def minimal_logger(name):
    return Logger(name)
