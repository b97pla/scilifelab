import shutil
import os
import logbook
import re
import yaml
import unittest

from ..classes import SciLifeTest
from classes import PmFullTest
from cement.core import handler
from scilifelab.pm.core.project import ProjectController, BcbioRunController
from scilifelab.utils.misc import filtered_walk

LOG = logbook.Logger(__name__)

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

j_doe_00_01 = os.path.abspath(os.path.join(filedir, "data", "production", "J.Doe_00_01"))
j_doe_00_04 = os.path.abspath(os.path.join(filedir, "data", "projects", "j_doe_00_04"))

SAMPLE="P001_102_index6"
SAMPLEFILE=os.path.join(filedir, "data", "projects", "j_doe_00_04", "samples.txt")

@unittest.skipIf(not os.getenv("MAILTO"), "not running production test: set $MAILTO environment variable to your mail address to test mailsend")
@unittest.skipIf(not os.getenv("DRMAA_LIBRARY_PATH"), "not running production test: no $DRMAA_LIBRARY_PATH")
class ProjectTest(PmFullTest):
    @classmethod
    def setUpClass(cls):
        if not os.getcwd() == filedir:
            os.chdir(filedir)
        LOG.info("Copy tree {} to {}".format(j_doe_00_01, j_doe_00_04))
        if not os.path.exists(j_doe_00_04):
            shutil.copytree(j_doe_00_01, j_doe_00_04)
        pattern = "-bcbb-config.yaml$"
        def yaml_filter(f):
            return re.search(pattern, f) != None
        yaml_files = filtered_walk(j_doe_00_04, yaml_filter)
        with open(SAMPLEFILE, "w") as fh:
            fh.write("\n".join(yaml_files[0:1]))

    def test_project_run(self):
        """Test running an amplicon project in project folder"""
        self.app = self.make_app(argv = ['project', 'run', 'j_doe_00_04', '--email', os.getenv("MAILTO"), '--restart', '--amplicon', '--drmaa', '-t', '00:05:00', '--partition', 'devel', '--sample', SAMPLE, '--debug', '--force', '--quiet'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProjectController)
        handler.register(BcbioRunController)
        self._run_app()
        os.chdir(filedir)

    def test_project_run_samplefile(self):
        """Test running project in project folder using a sample file as input."""
        self.app = self.make_app(argv = ['project', 'run', 'j_doe_00_04', '--restart', '--amplicon', '--sample', SAMPLEFILE, '--debug', '--force', '--quiet'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProjectController)
        handler.register(BcbioRunController)
        self._run_app()
        os.chdir(filedir)

    def test_project_run_only_failed(self):
        """Test running project in project folder on only failed samples. """
        psummary = os.path.join(j_doe_00_04, "P001_101_index3", "121015_BB002BBBXX", "project_summary.csv")
        if os.path.exists(psummary):
            os.unlink(psummary)
        self.app = self.make_app(argv = ['project', 'run', 'j_doe_00_04', '-n', '--only_failed', '--debug', '--force', '--quiet'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProjectController)
        handler.register(BcbioRunController)
        self._run_app()
        os.chdir(filedir)

    def test_project_monitor_job_fake_jobid(self):
        """Add fake JOBID and make sure command doesn't fail"""
        s = "P001_102_index6"
        jobidfile = os.path.join(j_doe_00_04, s, "121015_BB002BBBXX", "JOBID")
        LOG.info("Writing to jobid file {}".format(jobidfile))
        with open(jobidfile, "w") as fh:
            fh.write("1")
        self.app = self.make_app(argv = ['project', 'run', 'j_doe_00_04', '--drmaa', '--partition', 'devel', '-t', '00:10:00', '--sample', s, '--debug', '--force', '--quiet'], extensions=['scilifelab.pm.ext.ext_distributed'])
        handler.register(ProjectController)
        handler.register(BcbioRunController)
        self._run_app()
        os.chdir(filedir)
    
        
