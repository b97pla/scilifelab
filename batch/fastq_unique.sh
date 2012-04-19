#! /bin/bash

#SBATCH -A a2012043
#SBATCH -p devel
#SBATCH -t 1:00:00

#IN=$1
for IN in $*
do
  python ~/scilifelab/scripts/fastq_unique.py ${IN} &
done

