#!/bin/sh
# Returns Illumina Runs having Read3 "termination signal"
# Only works for PE reads/runs. A more precise detection scheme
# should be found on illumina_finished_msg.py on the BCBB pipeline

# This is just a quick commandline hint

find /mnt/illumina/*/ -maxdepth 1 -type f -iname 'Basecalling_Netcopy_complete_Read3.txt'
