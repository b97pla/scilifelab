"""
log module
"""
import sys
import logging
import bcbio.pipeline.config_loader as cl

def minimal_logger(namespace, debug=False):
    """Make and return a minimal console logger.

    NOTE: this does apparently *not* work with logbook as I first thought, and
    log handlers will *not* take care of output. If something is to be
    logged to a file in a module, the logger has to be implemented for
    that particular purpose.

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



def file_logger(namespace, config_file, log_file):
        CONFIG = cl.load_config(config_file)
        log_path = CONFIG['log_dir'] + '/' + log_file

        logger = logging.getLogger(namespace)
        logger.setLevel(logging.DEBUG)

        # file handler
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)

        # console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # formatter
        formatter = logging.Formatter("%(asctime)s (%(levelname)s) : %(message)s")
        fh.setFormatter(formatter)

        # add handlers to logger
        logger.addHandler(ch)
        logger.addHandler(fh)

        return logger
