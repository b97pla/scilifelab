#! /bin/bash

#SBATCH -A a2012043
#SBATCH -p core
#SBATCH -t 2:00:00

IN=$1
NAME=`basename ${IN}`
OUT=${NAME/_fastq.txt/-cutadapt.fastq}
LOG=${OUT/.fastq.gz/.cutadapt_metrics}

if [[ ${IN} =~ .*_1_fastq.* ]]
then
  # Adapter seq for read1 (reverse complement of common sequence of "NNN RNA PCR Primer Index N")
  #ADAPTER_SEQ=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCA
  # Adapter seq for read1 (32_Illumina RNA 3’ Adapter (RA3))
  ADAPTER_SEQ=TGGAATTCTCGGGTGCCAAGG
fi
if [[ ${IN} =~ .*_2_fastq.* ]]
then
  # Adapter seq for read1 (reverse complement of "100 Illumina RNA PCR Primer")
  ADAPTER_SEQ=TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT
fi

# If adapter sequence was passed, use that
if [[ ! -z $2 ]]
then
  ADAPTER_SEQ=$2 
fi

cutadapt -a ${ADAPTER_SEQ} -f fastq -O 5 -m 18 -o ${OUT} ${IN} >& ${LOG}

