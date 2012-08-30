"""
Test project subcontroller
"""
from cement.core import handler
from test_default import PmTest
from pmtools.controller.project import ProjectController

class InitProjectTest(PmTest):
    OUTPUT_FILES = []

    def setUp(self):
        super(InitProjectTest, self).setUp()
            
    def test_1_project_deliver(self):
        self.app = self.make_app(argv = ['project', 'deliver'])
        handler.register(ProjectController)
        self._run_app()
        
    def test_2_project_data_delivery(self):
        pass
