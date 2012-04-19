#! /bin/sh

#SBATCH -A a2010002
#SBATCH -p core
#SBATCH -t 2:00:00

python ~/scilifelab/scripts/split_demultiplexed.py $@
