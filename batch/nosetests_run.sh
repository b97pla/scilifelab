#!/bin/bash
#SBATCH -p devel
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J bcbb_testsuite
#SBATCH -A a2010002
#SBATCH -D /bubo/home/h5/roman/opt/bcbb/nextgen/tests
#SBATCH --mail-user=roman.valls.guimera@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -o nosetest_job.out
#SBATCH -e nosetest_job.err

nosetests -v -s
