"""
log module
"""
from logbook import Logger

log = Logger('scilifelab')
log.warn('This is a warning')

def minimal_logger(name, debug=False):
    return Logger(name)
