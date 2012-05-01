#! /bin/sh

#SBATCH -A a2010002
#SBATCH -p node
#SBATCH -t 2:00:00

REF=$1
P1=$2
P2=$3

t=8
if [ ! -z ${P2} ]
then
  t=4
fi

P1SAI=`basename ${P1}`
P1SAI=${P1SAI/.gz/.sai}

echo "Aligning ${P1}"
bwa aln -k 2 -l 18 -t ${t} -I ${REF} ${P1} > ${P1SAI} &

if [ ! -z ${P2} ]
then
  P2SAI=`basename ${P2}`
  P2SAI=${P2SAI/.gz/.sai}
  echo "Aligning ${P2}"
  bwa aln -k 2 -l 18 -t ${t} -I ${REF} ${P2} > ${P2SAI} &
fi

wait

SAM=${P1SAI/.sai/.sam}

echo "Generating alignments"
if [ -z ${P2} ]
then
  echo ""
  bwa samse ${REF} ${P1SAI} ${P1} > ${SAM}
else
  SAM=${SAM/_1-/-}
  bwa sampe ${REF} ${P1SAI} ${P2SAI} ${P1} ${P2} > ${SAM}
fi

echo "sam -> bam conversion"
BAM=${SAM/.sam/.bam}
samtools view -S -b ${SAM} > ${BAM}
SORTED=${BAM/.bam/.sort}
echo "Sorting bam"
samtools sort ${BAM} ${SORTED}
echo "Indexing bam"
samtools index ${SORTED}.bam

echo "Cleaning up"
rm ${P1SAI}
if [ ! -z ${P2} ]
then
  rm ${P2SAI}
fi
rm ${SAM}
rm ${BAM}


