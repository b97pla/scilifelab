#!/bin/bash
#SBATCH -p devel
#SBATCH -t 1:00:00
#SBATCH -J bcbb_testsuite
#SBATCH -A a2010002
#SBATCH --mail-user=pontus.larsson@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -o nosetest_job.out
#SBATCH -e nosetest_job.err

nosetests -v -s $1

