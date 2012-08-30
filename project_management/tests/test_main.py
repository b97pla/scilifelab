"""
Test main pm functionality
"""

from cement.utils import test
from test_default import PmTestApp

## NOTE: for some weird reason the help argument fails. Maybe
## because it is added as a default and thus not registered as an
## option?
class PmMainTest(test.CementTestCase):
    app_class = PmTestApp

    def test_1_config(self):
        app = self.make_app(argv=['--config'])
        try:
            app.setup()
            app.run()
        finally:
            app.close()

    ## Fails
    # def test_2_help(self):
    #     app = self.make_app(argv=['-h'])
    #     try:
    #         app.setup()
    #         app.run()
    #     finally:
    #         app.close()
