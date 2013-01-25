#! /bin/sh

#SBATCH -J fastqc
#SBATCH -A a2010002
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 2:00:00

F1=$1
RUNDIR=`dirname ${F1}`/..
OUTDIR=${RUNDIR}/fastqc
mkdir -p ${OUTDIR}
fastqc --outdir ${OUTDIR} $F1

