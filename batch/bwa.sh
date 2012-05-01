#! /bin/sh

set -e

FASTQ1=$1
FASTQ2=$2
REF=$3

REFFASTA=/bubo/nobackup/uppnex/reference/biodata/genomes/Hsapiens/hg19/bwa/hg19.fa

if [ ! -z ${REF} ]
then
  REFFASTA=${REF}
fi

###########################
## bwa aln ################
###########################

IDX1=`basename ${FASTQ1/.fastq.gz/-bwa-aln.sai}`
IDX2=`basename ${FASTQ2/.fastq.gz/-bwa-aln.sai}`

echo `date`" Running bwa aln" 
bwa aln -p 4 -q 20 ${REFFASTA} ${FASTQ1} > ${IDX1} &
bwa aln -p 4 -q 20 ${REFFASTA} ${FASTQ2} > ${IDX2}

wait

###########################
## bwa sampe ##############
###########################

SAMOUT=`basename ${FASTQ1/_1.fastq.gz/.sam}`

echo `date`" Running bwa sampe"
bwa sampe -r "${!RG}" -P ${REFDIR}/${REFFASTA} ${FASTQDIR}/${IDXFILE_1} ${FASTQDIR}/${IDXFILE_2} ${FASTQDIR}/${!FASTQFILE_1} ${FASTQDIR}/${!FASTQFILE_2} > ${BAMDIR}/${SAMOUT}

###########################
## sam to bam #############
###########################

BAMOUT=${SAMOUT/.sam/.bam}

echo "${SAMPLE}: Converting sam to bam"
samtools view -t ${REFDIR}/${REFFASTAIDX} -S -b ${BAMDIR}/${SAMOUT} > ${BAMDIR}/${BAMOUT}

###########################
## sort bam ###############
###########################

SORTOUT=${BAMOUT/.bam/-sort}

echo "${SAMPLE}: Sorting bam"
samtools sort ${BAMDIR}/${BAMOUT} ${BAMDIR}/${SORTOUT}
SORTOUT=${SORTOUT}.bam

###########################
## index bam ##############
###########################

echo "${SAMPLE}: Indexing bam"
samtools index ${BAMDIR}/${SORTOUT}

