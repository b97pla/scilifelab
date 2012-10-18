"""
log module
"""
import sys
from logbook import Logger
import logging


def minimal_logger(namespace, debug=False):
    """Make and return a minimal logger.

    NOTE: this does apparently *not* work with logbook as I first
    thought. Trick is to import logging and implement a logger that
    way. Log handlers will then take care of output.

    The current function is copied from cement.core.backend. 

    :param namespace: namspace of logger
    """

    log = logging.getLogger(namespace)
    formatter = logging.Formatter(
                "%(asctime)s (%(levelname)s) %(name)s : %(message)s")
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)   
    log.setLevel(logging.INFO)

    # FIX ME: really don't want to hard check sys.argv like this but
    # can't figure any better way get logging started (only for debug)
    # before the app logging is setup. Besides, this will fail for
    # tests since sys.argv will consist of the test call arguments.
    if '--debug' in sys.argv or debug:
        console.setLevel(logging.DEBUG)   
        log.setLevel(logging.DEBUG)

    log.addHandler(console)
    return log
