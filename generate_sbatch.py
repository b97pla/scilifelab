import sys
import os

if len(sys.argv) < 4:
    print "USAGE: python " + sys.argv[0] + " <path to folder with FASTQ files> <path to Bowtie index> <path to annotation GTF>"
    sys.exit(0)

old = False

fpath = sys.argv[1]
refpath = sys.argv[2]
annopath = sys.argv[3]
if len(sys.argv) == 5: 
    if sys.argv[4] == "-o"
    old = True
flist = os.listdir(fpath)

sample_names = []
read1forsample = {}
read2forsample = {}
for fname in flist:
    if 'fastq' not in fname: continue 
    read = fname.split("_")[-1]
    # print read
    tag = "_".join(fname.split("_")[3:-2])
    # 2_date_fcid_sample_1.fastq
    if not tag in sample_names: sample_names.append(tag)
    if read == "1.fastq": read1forsample[tag]=fname
    if read == "2.fastq": read2forsample[tag]=fname

print "Best guess for sample names: "
for n in sorted(sample_names):
    print n

r = raw_input("Press n to exit")
if r.upper() == "N": sys.exit(0)

for n in sorted(sample_names):
    print "Generating sbatch files for sample ", n
    oF = open("map_tophat_" + n + ".sh", "w")
    oF.write("#! /bin/bash -l\n")
    oF.write("#SBATCH -A a2010002\n")                   
    oF.write("#SBATCH -p node\n")
    oF.write("#SBATCH -t 35:00:00\n")
    oF.write("#SBATCH -J tophat_" + n + "\n")
    oF.write("#SBATCH -e tophat_" + n + ".err\n")
    oF.write("#SBATCH -o tophat_" + n + ".out\n")
    oF.write("#SBATCH --mail-user mikael.huss@scilifelab.se\n")
    oF.write("#SBATCH --mail-type=ALL\n")

    oF.write("module unload bioinfo-tools\n")
    oF.write("module unload samtools\n")
    oF.write("module unload tophat\n")
    oF.write("module unload cufflinks\n")
    oF.write("module unload htseq\n")

    oF.write("module load bioinfo-tools\n")
    oF.write("module load samtools/0.1.9\n")
    if old: oF.write("module load tophat/1.0.14\n")
    else: oF.write("module load tophat/1.3.3\n")
    oF.write("module load cufflinks/1.2.1\n")
    oF.write("module load htseq/0.5.1\n")

    # TopHat

    size = raw_input("Fragment size including adapters for sample" + n + "? ")
    innerdist = int(size) - 320
    if not size.isdigit(): sys.exit(0)

    oF.write("tophat -o tophat_out_" + n + " --solexa1.3-quals -p 8 -r " + str(innerdist) + " " + refpath + " " + fpath + "/" + read1forsample[n] + " " + fpath + "/" + read2forsample[n] + "\n")
    # Samtools -> Picard

    oF.write("cd tophat_out_" + n + "\n")
    if old: oF.write("samtools view -bT concat.fa.fa -o accepted_hits_" + n + ".bam " + "accepted_hits.sam\n")
    else: oF.write("mv accepted_hits.bam accepted_hits_" + n + ".bam\n")
    oF.write("java -Xmx2g -jar /home/lilia/glob/src/picard-tools-1.29/SortSam.jar INPUT=accepted_hits_" + n + ".bam OUTPUT=accepted_hits_sorted_" + n + ".bam SORT\
_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT\n")
    oF.write("java -Xmx2g -jar /home/lilia/glob/src/picard-tools-1.29/MarkDuplicates.jar INPUT=accepted_hits_sorted_" + n + ".bam OUTPUT=accepted_hits_sorted_dup\
Removed_" + n + ".bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=" + n + "_picardDup_metrics VALIDATION_STRINGENCY=LENIENT\n")

    # HTSeq
    oF.write("samtools view accepted_hits_sorted_dupRemoved_" + n + ".bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_" + n + ".sam\n")
    oF.write("python -m HTSeq.scripts.count -s no -q accepted_hits_sorted_dupRemoved_prehtseq_" + n + ".sam " + annopath + " > " + n + ".counts\n")

    # Cufflinks

    oF.write("samtools view accepted_hits_sorted_dupRemoved_" + n + ".bam |sort -k 3,3 -k 4,4n > accepted_hits_sorted_dupRemoved_col34Sorted_" + n + ".sam\n")
    oF.write("cufflinks -p 8 -G " + annopath + " -o cufflinks_out_" + n + " accepted_hits_sorted_dupRemoved_col34Sorted_" + n + ".sam\n")

    # Clean up
    oF.write("rm accepted_hits_sorted_" + n + ".bam\n")
    oF.write("rm accepted_hits_sorted_dupRemoved_prehtseq_" + n + ".sam\n")
    oF.write("rm accepted_hits_sorted_dupRemoved_col34Sorted_" + n + ".sam\n")

    oF.close()
