#! /bin/sh

#SBATCH -A a2010002
#SBATCH -p core
#SBATCH -t 1:00:00

md5sum $1 >> $2

