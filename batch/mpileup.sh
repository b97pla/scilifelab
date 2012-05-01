#! /bin/sh

#SBATCH -A a2010002
#SBATCH -p core
#SBATCH -t 24:00:00
#SBATCH -J mpileup

REF=/bubo/nobackup/uppnex/reference/biodata/genomes/Hsapiens/hg19/seq/hg19.fa

for IN in $@
do
  RAWOUT=${IN/.bam/-mpileup-raw.bcf}
  FILTEREDOUT=${RAWOUT/raw.b/filtered.v}
  PILEOUT=${RAWOUT/.bcf/.pileup}

  samtools mpileup -uf ${REF} ${IN} | bcftools view -bvcg - > ${RAWOUT}  
  bcftools view ${RAWOUT} | vcfutils.pl varFilter -D100 -d 10 > ${FILTEREDOUT}
  rm ${RAWOUT}  
done

