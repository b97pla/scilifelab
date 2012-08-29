"""
Test main pm functionality
"""

from cement.utils import test
from pmtools import PmController, PmApp
from test_default import config_defaults, PmTestApp

## NOTE: for some weird reason the help argument fails. Maybe
## because it is added as a default and thus not registered as an
## option?
class PmMainTest(test.CementTestCase):
    app_class = PmApp

    def test_1_config(self):
        app = PmApp('pm', base_controller=PmController, config_defaults = config_defaults, extensions = ['json'], argv=['--config'], config_files=[])
        try:
            app.setup()
            app.run()
        finally:
            app.close()
