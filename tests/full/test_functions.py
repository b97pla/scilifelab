import os
import yaml
import logbook
from scilifelab.bcbio import prune_pp_platform_args
from scilifelab.experiment.project import find_samples

from ..classes import SciLifeTest

LOG = logbook.Logger(__name__)

filedir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

j_doe_00_01 = os.path.abspath(os.path.join(filedir, "data", "production", "J.Doe_00_01"))
SAMPLE = "P001_101_index3"
ALLSAMPLES = [SAMPLE, "P001_102_index6"]
FLOWCELL = '120924_AC003CCCXX'

class BcbioFunctionTest(SciLifeTest):
    """Test functions in scilifelab.bcbio"""
    def test_prune_pp_platform_args(self):
        ppfile = os.path.join(j_doe_00_01, SAMPLE, FLOWCELL, "{}-post_process.yaml".format(SAMPLE))
        with open(ppfile, "r") as fh:
            conf = yaml.load(fh)
        newconf = prune_pp_platform_args(conf, keep_opts=["-t"])
        self.assertIn("-t", newconf['distributed']['platform_args'].split())
        self.assertNotIn("-A", newconf['distributed']['platform_args'].split())


class ExperimentFunctionTest(SciLifeTest):
    """Test functions in scilifelab.experiment"""
    def test_find_samples(self):
        plist = find_samples(j_doe_00_01)
        self.assertEqual([os.path.basename(x) for x in plist], ALLSAMPLES)
        plist = find_samples(j_doe_00_01, sample="P111")
        self.assertEqual([os.path.basename(x) for x in plist], [])
        plist = find_samples(j_doe_00_01, sample=SAMPLE)
        self.assertEqual([os.path.basename(x) for x in plist], [SAMPLE])
        plist = find_samples(j_doe_00_01, sample="P001")
        self.assertEqual([os.path.basename(x) for x in plist], [])
        samplefile = os.path.join(j_doe_00_01, "samples.txt")
        with open(samplefile, "w") as fh:
            fh.write(SAMPLE)
        plist = find_samples(j_doe_00_01, sample=samplefile)
        self.assertEqual([os.path.basename(x) for x in plist], [SAMPLE])
