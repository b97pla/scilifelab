#! /bin/bash -l

#SBATCH -A b2013064
#SBATCH -o mapTopHat.out
#SBATCH -e mapTopHat.err
#SBATCH -J mapTopHat.job 
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH --mail-user pelin.akan@scilifelab.se
#SBATCH --mail-type=ALL


# Arguments:  sbatch mapTopHat.sh samplefolder
#
# Command:  tophat -p 8  --library-type fr-firststrand --no-coverage.search -o samplename reference_index forwardreads reversereads
#
# Result: sorted BAM file   and some statistics

module load bioinfo-tools
module load tophat/2.0.4
module load picard
module load rseqc

#This needs to change every time
cd /bubo/home/h20/pelin/app_b2013064/nobackup/pelin/projects/M.Jagodic_12_01_2ndround/analysis
#do not forget to soft link all the folders of the samples in the following folder
#do not forget to soft link the genome and use the right genome for the mappings!
 
#Name of the sample is the argument, each sample has its own folder (which needs to be soft-linked) under that folder, each sample lane has a seperate folder

samplename=$1
echo $samplename  
#Generate an array containing all the fastq files belonging to the sample
space=","; 
for reads in `ls $samplename/*/*_1.*`; do  readlist1+=($reads$space); done;
for reads in `ls $samplename/*/*_2.*`; do  readlist2+=($reads$space); done;

#get the array size
asize1=${#readlist1[@]}
echo $asize1
asize2=${#readlist2[@]}
echo $asize2

#if there is more than one forward or reverse reads, enter them comma-seperated
#get rid of the last comma
lastread1=`echo ${readlist1[$asize1-1]}| cut -d',' -f 1`
lastread2=`echo ${readlist2[$asize2-1]}| cut -d',' -f 1`
unset readlist1[$asize1-1]
unset readlist2[$asize2-1]
readlist1=("${readlist1[@]}" $lastread1)
readlist2=("${readlist2[@]}" $lastread2)

#get rid of the space between commas
for i in ${readlist1[@]};do forwardreads+=$i; done;
for j in ${readlist2[@]};do reversereads+=$j; done;


echo "Map with Tophat"
echo tophat -p 8 --library-type fr-firststrand --no-coverage-search -o $samplename rn4 $forwardreads $reversereads
tophat -p 8 --library-type fr-firststrand --no-coverage-search -o $samplename rn4 $forwardreads $reversereads

echo "Mark Duplicates"
java -Xmx2g -jar /bubo/sw/apps/bioinfo/picard/1.92/kalkyl/MarkDuplicates.jar I=$samplename/accepted_hits.bam O=$samplename\.sorted.dup.bam METRICS_FILE=$samplename\.dupmetrics.txt REMOVE_DUPLICATES=0 ASSUME_SORTED=1 VALIDATION_STRINGENCY=SILENT

echo "Check Strand Specificity"
infer_experiment.py -i $samplename\.sorted.dup.bam -r ../../../annotations/rn4_RGDGenes.bed >$samplename\.strandspec.txt

echo "Generate Complexity Curve"
lc_extrap -v -B $samplename\.sorted.dup.bam -o $samplename\.ccurve.txt

echo "Run BAMQC"
qualimap --java-mem-size=1000M bamqc -bam $samplename\.sorted.dup.bam -c -outdir $samplename\_QC/


