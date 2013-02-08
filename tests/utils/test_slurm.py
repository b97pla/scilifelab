"""Test the utils/slurm.py functionality
"""
import subprocess
import unittest
from mock import Mock

import scilifelab.utils.slurm as sq
from scilifelab.pm.ext.ext_distributed import convert_to_drmaa_time

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
        
        

class TestDrmaa(unittest.TestCase):

    def test_drmaa_time_string(self):
        """Test parsing of time string formatted as d-hh:mm:ss and translate days to hours"""
        t_new = convert_to_drmaa_time("4-10:00:00")
        self.assertEqual(t_new, "106:00:00")
        t_new = convert_to_drmaa_time("10:00:00")
        self.assertEqual(t_new, "10:00:00")
        t_new = convert_to_drmaa_time("3:00:00")
        self.assertEqual(t_new, "03:00:00")
        t_new = convert_to_drmaa_time("10:00")
        self.assertEqual(t_new, "00:10:00")
        t_new = convert_to_drmaa_time("0:00")
        self.assertEqual(t_new, "00:00:00")
        t_new = convert_to_drmaa_time("144:00:00")
        self.assertEqual(t_new, "144:00:00")
