#!/bin/sh
## Script ready to be put on python virtualenv environment
## should be called from ~/.virtualenv/<workspace>/bin/postactivate
## *NOT* from .bashrc or 

# No python modules are loaded since only the ones
# in the virtual environment will be used.
# "mkvirtualenv --no-site-packages <workspace>" is assumed

MODS="bioinfo-tools
bowtie/0.12.7
samtools
bwa/0.5.9
maq
mosaik-aligner
python/2.6.6
CASAVA
picard
GATK
cufflinks
tophat
Fastx
FastQC
R/2.12.2
emboss
gnuplot
BEDTools"

module purge

for mod in $MODS
do
	module load $mod >& /dev/null
done

# We unset PYTHONHOME set by the module system, 
# otherwise the system will not use the python
# of virtualenv
unset PYTHONHOME
