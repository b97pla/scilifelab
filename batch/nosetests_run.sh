#!/bin/bash
#SBATCH -p devel
#SBATCH -N 1
#SBATCH -t 00:60:00
#SBATCH -J bcbb_testsuite
#SBATCH -A a2010002
#SBATCH -D .
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o nosetest_job.out
#SBATCH -e nosetest_job.err

nosetests -v -s $1
