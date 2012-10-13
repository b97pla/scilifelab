"""
log module
"""
from logbook import Logger

def minimal_logger(namespace):
    """Make and return a minimal logger

    :param namespace: namspace of logger
    """
    return Logger(namespace)
