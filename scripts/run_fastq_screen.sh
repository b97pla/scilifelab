#! /bin/sh

#SBATCH -J fq_screen
#SBATCH -A a2010002
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 3:00:00

F1=$1
F2=${F1/1_fastq.txt/2_fastq.txt}
RUNDIR=`dirname ${F1}`/..
#RUNDIR=${PWD}
OUTDIR=${RUNDIR}/fastq_screen
mkdir -p ${OUTDIR}

#fastq_screen --outdir ${RUNDIR}/fastq_screen --subset 2000000 --illumina --multilib $F1
fastq_screen --outdir ${OUTDIR} --subset 2000000 --illumina --multilib $F1 --paired $F2

