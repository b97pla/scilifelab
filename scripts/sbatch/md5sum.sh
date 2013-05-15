#! /bin/sh

#SBATCH -A a2010002
#SBATCH -p core
#SBATCH -t 6:00:00
#SBATCH -J md5

md5sum $1 >> $2

