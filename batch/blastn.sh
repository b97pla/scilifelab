#! /bin/sh

#SBATCH -p devel
#SBATCH -A a2010002
#SBATCH -t 1:00:00
#SBATCH -J blastn

DB=$1
QFILE=$2

FNAME=`basename ${QFILE}`
OUTFILE=${FNAME/.fastq.gz/.blastn.out}

zcat ${QFILE} | fastq_to_fasta -n | blastall -p blastn -m 8 -F F -a 8 -W 18 -q 5 -r -4 -G 6 -E 5 -d ${DB} -o ${OUTFILE}
 