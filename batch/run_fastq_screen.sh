#! /bin/sh

#SBATCH -J fq_screen
#SBATCH -A a2010002
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH --qos=seqver

F1=$1
F2=${F1/1_fastq.txt/2_fastq.txt}
RUNDIR=`dirname ${F1}`/..
OUTDIR=${RUNDIR}/fastq_screen
mkdir -p ${OUTDIR}

if [ ${F1##*.} == "gz" ]
then
  F1=${F1/.gz/}
  gunzip -c ${F1}.gz > ${F1}
fi

if [ ${F2##*.} == "gz" ]
then
  F2=${F2/.gz/}
  gunzip -c ${F2}.gz > ${F2}
fi

#fastq_screen --outdir ${OUTDIR} --subset 2000000 --illumina --multilib $F1
fastq_screen --outdir ${OUTDIR} --subset 2000000 --illumina --multilib $F1 --paired $F2

