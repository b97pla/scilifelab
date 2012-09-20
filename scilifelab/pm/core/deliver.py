"""
Pm deliver module

Perform data delivery.

Commands:
	miseq	deliver miseq data from archive to customer
        
"""

import os
import shutil

from cement.core import controller
from scilifelab.pm.core.controller import AbstractBaseController

## Main delivery controller
class DeliveryController(AbstractBaseController):
    """
    Functionality for deliveries
    """
    class Meta:
        label = 'deliver'
        description = 'Deliver data'
        arguments = [
            (['runid'], dict(help="run id")),
            (['projectid'], dict(help="project id")),
            ]

    @controller.expose(hide=True)
    def default(self):
     	print __doc__

    @controller.expose(help="Deliver miseq data from RUNID to PROJECTID")
    def miseq(self):
        self._not_implemented("On hold. Archive data needs to be processed first")
        assert os.path.exists(os.path.join(self.config.get("archive", "root")), self.pargs.runid), "no such runid " % self.pargs.runid
        assert os.path.exists(os.path.join(self.config.get("archive", "root")), self.pargs.runid), "no such project id" % self.pargs.projectid

    #@controller.expose(help="Custom delivery. Deliver ")
