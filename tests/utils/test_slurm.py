"""Test the utils/slurm.py functionality
"""
import subprocess
import unittest
from mock import Mock

import scilifelab.utils.slurm as sq

class TestSlurm(unittest.TestCase):
    
    
    def test__get_slurm_jobid(self):
        """Extract the jobid for a slurm job name
        """
        
        # Mock the system calls
        subprocess.check_output = Mock(return_value='')
        # Assert that getting non-existing jobs return an empty job list
        self.assertListEqual([],sq.get_slurm_jobid("jobname"),
                             "Querying for jobid of non-existing job should return an empty list")
        # Assert that a returned job id is parsed correctly
        for jobids in [[123456789],[123456789,987654321]]:
            subprocess.check_output = Mock(return_value="\n".join([str(jid) for jid in jobids]))
            self.assertListEqual(jobids,sq.get_slurm_jobid("jobname"),
                                 "Querying for jobid of existing job did not return the correct value")
        
        