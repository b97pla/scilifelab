#!/usr/bin/env python
import sys
import os
import optparse
import glob
import yaml
import couchdb
import bcbio.google
import bcbio.pipeline.config_loader as cl

usage = """
Generate sbatch files for running an mRNA-seq pipeline.
Usage:

generate_sbatch.py <Flow cell ID, eg 121113_BD1HG4ACXX> <project ID> <Bowtie index for reference genome> <gene/transcript annotation GTF file>
The four first options are mandatory. There are further flags you can use:

-c, --config: Provide config file (post_process.yaml)
-o, --old: Run old TopHat (1.0.14) instead - mainly to reproduce old runs
-t, --projtag: Provide a project tag that will be shown in the queuing system to distinguish from other TopHat runs
-p, --phred64: Use phred64 quality scale
-f, --fai: Generate UCSC Genome Browser compatible BigWig tracks, using the specified FASTA index (fai) file
-m, --mail: Specify a mailing address for SLURM
-a, --alloc-time: Time to allocate in SLURM. Hours:minutes:seconds or days-hours:minutes:seconds
-w, --fpath: Path to fastq files. If not given, the script assumes the standars structure of the analysis directory. And that you run the script from the intermediate directory.

"""

if len(sys.argv) < 5:
    sys.exit(usage)

##	Couchdb functions
def find_proj_from_view(proj_db, proj_id):
        view = proj_db.view('project/project_name')
        for proj in view:
		if proj.key == proj_id:
                	return proj.value
	return None

def get_size(n,info):
	preps = 'BCDE'
	prep = ''
	name = n.replace('-', '_').replace(' ', '').split("_index")[0].split("_ss")[0].split("_dual")[0].strip().strip('F')
        while name[-1] in preps:
        	prep = name[-1] + prep
               	name = name[0: -1]
	if prep == '':
		prep = 'A'
	for samp in info['samples']:
        	if samp==name:
                	size=info['samples'][samp]['library_prep'][prep]['average_size_bp']
		elif name[0]=='P' and name.split('_')[1]==samp:
			print "Best gues for sample found on statusDB corresponding to sample name ",n ," is ",samp
			r = raw_input("If this seems to be correct press y	")
			if r.upper() == "Y":
				size=info['samples'][samp]['library_prep'][prep]['average_size_bp']
	return int(float(size))


parser = optparse.OptionParser()
parser.add_option('-c', '--config', action="store", dest="conffile", default=os.path.expanduser("~/opt/config/post_process.yaml"), help="Specify config file (post_process.yaml)")
parser.add_option('-o', '--old', action="store_true", dest="old", default="False", help="Run old TopHat (1.0.14) instead - mainly to reproduce old runs")
parser.add_option('-t', '--projtag', action="store", dest="projtag", default="", help="Provide a project tag that will be shown in the queuing system to distinguish from other TopHat runs")
parser.add_option('-p', '--phred64', action="store_true", dest="phred64", default="False", help="Use phred64 quality scale")
parser.add_option('-f', '--fai', action="store", dest="fai", default="", help="Provide FASTA index file for generating UCSC bigwig tracks")
parser.add_option('-m', '--mail', action="store", dest="mail", default="mikael.huss@scilifelab.se", help="Specify a mailing address for SLURM mail notifications")
parser.add_option('-a', '--alloc-time', action="store", dest="hours", default="40:00:00", help="Time to allocate in SLURM. Please specify as hours:minutes:seconds or days-hours:minutes:seconds")
parser.add_option('-w', '--fpath', dest="fpath", default=False, help="Path to fastq files. If not given, the script assumes the standars structure of the analysis directory'")
(opts, args) 	= parser.parse_args()
old		= opts.old
phred64 	= opts.phred64
fai 		= opts.fai
projtag 	= opts.projtag
mail 		= opts.mail
hours 		= opts.hours
conffile	= opts.conffile
fpath 		= opts.fpath


proj_ID         = sys.argv[2]
refpath         = sys.argv[3]
annopath        = sys.argv[4]

conf            = cl.load_config(conffile)

quantifyer_version      = conf['custom_algorithms']['RNA-seq analysis']['quantifyer_version']
counts_version          = conf['custom_algorithms']['RNA-seq analysis']['counts_version']
dup_remover_version     = conf['custom_algorithms']['RNA-seq analysis']['dup_remover_version']
picard_version		= conf['custom_algorithms']['RNA-seq analysis']['picard_version']
picard_tools		= conf['custom_algorithms']['RNA-seq analysis']['picard_tools']
BEDTools_version	= conf['custom_algorithms']['RNA-seq analysis']['BEDTools_version']
bedGraphToBigWig	= conf['custom_algorithms']['RNA-seq analysis']['bedGraphToBigWig']
URL             	= conf['couch_db']['maggie_url']
couch           	= couchdb.Server("http://" + URL)
proj_db         	= couch['projects']
key             	= find_proj_from_view(proj_db, proj_ID)
info            	= proj_db[key]


if not len ( hours.split(':') ) == 3:
	sys.exit("Please specify the time allocation string as hours:minutes:seconds or days-hours:minutes:seconds")

if phred64 == True:
	qscale = '--solexa1.3-quals'
else:
	qscale = ''

if not fpath:
	p = os.getcwd()
	fpath = p.split('intermediate')[0] + 'data/' + sys.argv[1]
print fpath

test=flist[0].split('.')[-1]
if not (test=='fastq')|(test=='gz'):
	sys.exit('Something wrong with the path to the fastqfiles: '+ fpath)

refpath = sys.argv[3]
annopath = sys.argv[4]

sample_names = []
read1forsample = {}
read2forsample = {}

try:
    conf = yaml.load(open(conffile))
except:
    sys.exit("Could not open configuration file " + conffile)

for fname in flist:
    read = fname.split("_")[-1]
    tag = "_".join(fname.split("/")[-1].split("_")[3:-2])
    print fname.split("_")
    # 2_date_fcid_sample_1.fastq
    if not tag in sample_names: sample_names.append(tag)
    if (read == "1.Q25.fastq") | (read == "1.fastq") | (read == "1.fastq.gz") | (read == "1.Q25.fastq.gz"): read1forsample[tag]=fname.split('/')[-1]
    if (read == "2.Q25.fastq") | (read == "2.fastq") | (read == "2.fastq.gz") | (read == "2.Q25.fastq.gz"): read2forsample[tag]=fname.split('/')[-1]
    print fname.split('/')[-1]
print "Best guess for sample names: "
for n in sorted(sample_names):
    print n
print read2forsample.keys()
r = raw_input("Press n to exit")
if r.upper() == "N": sys.exit(0)

#sFile = open("sample_names.txt", "w")
#for n in sorted(sample_names):
#    sFile.write(n + "\n")

for n in sorted(sample_names):
    print "Generating sbatch files for sample ",n
    oF = open("map_tophat_" + n + ".sh", "w")
    oF.write("#! /bin/bash -l\n")
    oF.write("#SBATCH -A a2012043\n")
    oF.write("#SBATCH -p node\n")
    oF.write("#SBATCH -t " + hours + "\n")
    oF.write("#SBATCH -J tophat_" + n + projtag + "\n")
    oF.write("#SBATCH -e tophat_" + n + projtag + ".err\n")
    oF.write("#SBATCH -o tophat_" + n + projtag + ".out\n")
    oF.write("#SBATCH --mail-user=" + mail + "\n")
    oF.write("#SBATCH --mail-type=ALL\n")
    oF.write("module unload bioinfo-tools\n")
    #oF.write("module unload samtools\n")
    oF.write("module unload tophat\n")
    oF.write("module unload cufflinks\n")
    oF.write("module unload htseq\n")

    oF.write("module load bioinfo-tools\n")
    # oF.write("module load samtools/0.1.9\n")

    tool_versions = conf['custom_algorithms']['RNA-seq analysis']

    if old == True: oF.write("module load tophat/1.0.14\n")
    else: oF.write("module load tophat/" + tool_versions['aligner_version'] + "\n")
    oF.write("module load cufflinks/" + tool_versions['quantifyer_version'] + "\n")
    oF.write("module load htseq/" + tool_versions['counts_version'] + "\n")

    # TopHat
    try:
	# trying to get average fragemnt length from couchDB
	size=get_size(n,info)
	print "Average fragment length ",str(size)
    except:
	print ""
    	size = raw_input("Could not find information on statusDB. Enter manualy fragment size including adapters for sample " + n + ": ")
	pass
    innerdist = int(size) - 101 - 101 - 121
    size=str(size)
    if not size.isdigit(): sys.exit(0)
    print fpath
    oF.write("tophat -o tophat_out_" + n + " " + qscale + " -p 8 -r " + str(innerdist) + " " + refpath + " " + fpath + "/" + read1forsample[n] + " " + fpath + "/" + read2forsample[n] + "\n")
    # Samtools -> Picard

    oF.write("cd tophat_out_" + n + "\n")
    if old == True: oF.write("samtools view -bT concat.fa.fa -o accepted_hits_" + n + ".bam " + "accepted_hits.sam\n")
    else: oF.write("mv accepted_hits.bam accepted_hits_" + n + ".bam\n")
    oF.write("java -Xmx2g -jar /home/lilia/glob/src/picard-tools-" + tool_versions['dup_remover_version'] + "/SortSam.jar INPUT=accepted_hits_" + n + ".bam OUTPUT=accepted_hits_sorted_" + n + ".bam SORT\
_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT\n")
    oF.write("java -Xmx2g -jar /home/lilia/glob/src/picard-tools-" + tool_versions['dup_remover_version'] + "/MarkDuplicates.jar INPUT=accepted_hits_sorted_" + n + ".bam OUTPUT=accepted_hits_sorted_dup\
Removed_" + n + ".bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=" + n + "_picardDup_metrics VALIDATION_STRINGENCY=LENIENT\n")

    # HTSeq
    oF.write("samtools view accepted_hits_sorted_dupRemoved_" + n + ".bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_" + n + ".sam\n")
    oF.write("python -m HTSeq.scripts.count -s no -q accepted_hits_sorted_dupRemoved_prehtseq_" + n + ".sam " + annopath + " > " + n + ".counts\n")

    # Cufflinks

    # To use indexed BAM file in Cufflinks

    oF.write("samtools index accepted_hits_sorted_dupRemoved_" + n + ".bam\n")
    oF.write("cufflinks -p 8 -G " + annopath + " -o cufflinks_out_" + n + " accepted_hits_sorted_dupRemoved_" + n + ".bam\n")

    # To use col 3-4 sorted SAM file in Cufflinks

p=os.getcwd()

if not os.path.exists(fpath):
	sys.exit('no such dir '+fpath)


flist = glob.glob(str(fpath + '/*'))
file_info={}

for f in flist:
	fname	=	f.split('/')[-1]
	if (fname.split('.')[-1] == "fastq") | ( fname.split('.')[-2] == "fastq"):
    		lane_run 	= "_".join(fname.split("_")[0:3])
    		tag		= "_".join(fname.split("_")[3:-2])
    		if not file_info.has_key(lane_run):
    			file_info[lane_run]={}
    		if not file_info[lane_run].has_key(tag):
			file_info[lane_run][tag]=[]
		file_info[lane_run][tag].append(fname)

print file_info.keys()

print "Best guess for sample names: "
for lane in file_info:
	print lane
	print sorted(file_info[lane].keys())

r = raw_input("Press n to exit")
if r.upper() == "N": sys.exit(0)



for lane in file_info:
    an_path = p+'/'+lane
    if not os.path.exists(an_path):
	os.mkdir(an_path)

    print an_path
    for samp in sorted(file_info[lane]):
	## Handeling some anoying options in a ugly way
	if old == True:
        	aligner_version = "1.0.14"
        	sam_bam="samtools view -bT concat.fa.fa -o accepted_hits_"+samp+".bam accepted_hits.sam"
	else:
        	aligner_version = conf['custom_algorithms']['RNA-seq analysis']['aligner_version']
        	sam_bam="mv accepted_hits.bam accepted_hits_"+samp+".bam"
	if fai != '':
        	make_fai="""genomeCoverageBed -bga -split -ibam accepted_hits_sorted_dupRemoved_"""+samp+""".bam -g """+fai+""" > sample_"""+samp+""".bga
"""+ bedGraphToBigWig +""" sample_"""+samp+""".bga"""+fai+""" sample_"""+samp+""".bw"""
	else:
		make_fai=''

	## trying to get average fragemnt length from couchDB
        try:
                size = get_size(samp,info)
                print "Average fragment length ", str(size)
        except:
                print ""
                size = raw_input("Could not find information on statusDB. Enter manualy fragment size including adapters for sample " + samp + ": ")
                pass
        innerdist = int(size) - 101 - 101 - 121
        size=str(size)
        if not size.isdigit(): sys.exit(0)

	### Generating sbatch file for sample
	f = open(an_path+"/map_tophat_"+samp+".sh", "w")
	print >>f, """#! /bin/bash -l

#SBATCH -A a2012034
#SBATCH -p node
#SBATCH -t {0}
#SBATCH -J tophat_{1}
#SBATCH -e tophat_{1}.err
#SBATCH -o tophat_{1}.out
#SBATCH --mail-user={2}
#SBATCH --mail-type=ALL

module unload bioinfo-tools
module unload bowtie
module unload tophat
module unload cufflinks
module unload htseq
module unload picard

module load bioinfo-tools
module load samtools/0.1.9
module load tophat/{3}
module load cufflinks/{4}
module load htseq/{5}
module load picard/{16}
module load BEDTools/{18}

tophat -o tophat_out_{1} {15} -p 8 -r {6} {7} {8}/{9} {8}/{10}
cd tophat_out_{1}
{13}
java -Xmx2g -jar {17}/SortSam.jar INPUT=accepted_hits_{1}.bam OUTPUT=accepted_hits_sorted_{1}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
java -Xmx2g -jar {17}/MarkDuplicates.jar INPUT=accepted_hits_sorted_{1}.bam OUTPUT=accepted_hits_sorted_dupRemoved_{1}.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE={1}_picardDup_metrics VALIDATION_STRINGENCY=LENIENT
samtools view accepted_hits_sorted_dupRemoved_{1}.bam |sort > accepted_hits_sorted_dupRemoved_prehtseq_{1}.sam
python -m HTSeq.scripts.count -s no -q accepted_hits_sorted_dupRemoved_prehtseq_{1}.sam {12}  > {1}.counts
samtools index accepted_hits_sorted_dupRemoved_{1}.bam
cufflinks -p 8 -G {12} -o cufflinks_out_{1} accepted_hits_sorted_dupRemoved_{1}.bam

{14}

rm accepted_hits_sorted_{1}.bam
rm accepted_hits_sorted_dupRemoved_prehtseq_{1}.sam""".format(hours ,samp ,mail ,aligner_version,quantifyer_version,counts_version,str(innerdist),refpath,fpath,file_info[lane][samp][0],file_info[lane][samp][1],dup_remover_version,annopath,sam_bam,make_fai, qscale,picard_version,picard_tools,BEDTools_version)
    	f.close()

