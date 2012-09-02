"""dry"""

import sys

def dry(ctrl, message, func, *args, **kw):
    if ctrl.pargs.dry_run:
        print >> sys.stderr, "(DRY_RUN): " + message
        ctrl.app._output_data["stderr"].write("(DRY_RUN): " + message)
        return
    ctrl.log.info(message)
    return func(*args, **kw)
