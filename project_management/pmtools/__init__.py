"""
Pipeline Management Tools

Usage: pm command [options]
"""

from cement.core import foundation, controller, handler

## Abstract base controller -- for sharing arguments
class AbstractBaseController(controller.CementBaseController):
    class Meta:
        arguments = [
            (['-n', '--dry_run'], dict(help="dry_run - don't actually do anything")),
            ]

    def _setup(self, base_app):
        super(AbstractBaseController, self)._setup(base_app)
        self.shared_config = dict()


## Main pm base controller
class PmController(controller.CementBaseController):
    class Meta:
        label = 'pm'
        description = ''

    @controller.expose(hide=True)
    def default(self):
        print dir(self)





