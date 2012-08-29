"""
Test main pm functionality
"""

from cement.utils import test
from pmtools import PmController, PmApp
from test_default import config_defaults, PmTestApp

class PmMainTest(test.CementTestCase):
    app_class = PmTestApp

    def test_1_help(self):
        app = self.make_app(argv=['--help'])
        app.setup()
        print dir(app)
        print app.argv
        print dir(app.config)
        print dir(app._meta)
        print app._meta.label
        print app._meta.base_controller.Meta.label
        app.run()
        app.close()

