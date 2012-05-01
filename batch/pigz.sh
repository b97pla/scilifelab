#! /bin/sh

#SBATCH -A a2012043
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J pigz

IN=$1
NAME=`basename ${IN}`
OUT=${NAME}.gz
cat ${IN} | pigz -p 1 -c > ${OUT}
