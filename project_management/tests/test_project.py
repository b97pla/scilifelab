"""
Test project subcontroller
"""
from cement.core import handler
from cement.utils import test
from test_default import PmTestApp
from pmtools.controller.project import ProjectController

class InitProjectTest(test.CementTestCase):
    app_class = PmTestApp

    def setUp(self):
        super(InitProjectTest, self).setUp()
            
    def test_1_project_deliver(self):
        app = self.make_app(argv = ['project', 'deliver'])
        handler.register(ProjectController)
        app.setup()
        app.run()
        app.close()
        
    def test_2_project_data_delivery(self):
        pass
