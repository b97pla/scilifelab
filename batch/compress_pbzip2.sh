#!/bin/bash
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -J pbzip2
#SBATCH -A a2010002
#SBATCH --mail-user=per.unneberg@scilifelab.se
#SBATCH --mail-type=ALL

# Compress everything:
# http://www.dataspora.com/2010/08/the-seven-secrets-of-successful-data-scientists/
# ... using all cores with pbzip2

function usage {
    echo "Usage: $0 <-d> <-s> [dir|file] <dest_dir>"
    echo "Example: sbatch compress_slurm.sh bigfolder /store/longterm/mount/"
    echo "Output: file.bz2 if it's a single file or dir.tar.bz2, -d for decompressing"
}


if [ $# == 0 ]; then
	usage && exit
fi

# Decompress
if [ $1 == "-d" ]; then
	tar xf $2.tar.bz2 --use-compress-prog=pbzip2
	exit
fi

# If its a directory, tar it too
if [ -d $1 ]
then
    tar cf $2`basename $1`.tar.bz2 --use-compress-prog=pbzip2 $1
elif [ -f $1 ]
then
# In this case, it bypasses core autodetection and sets 16 cores
# (16 "virtual cores" by HyperThreading):
    pbzip2 -p16 $1
else
    usage
fi
