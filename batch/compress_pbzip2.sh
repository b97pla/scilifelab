#!/bin/bash
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5:00:00
#SBATCH -J pbzip2
#SBATCH -A YOUR_PROJECT_OR_QUEUE
#SBATCH --mail-user=YOURMAIL@example.com
#SBATCH --mail-type=ALL

# Compress everything:
# http://www.dataspora.com/2010/08/the-seven-secrets-of-successful-data-scientists/
# ... using all cores with pbzip2

function usage {
    echo "Usage: compress_pbzip2.sh [dir|file]"
    echo "Output: file.bz2 if it's a single file or dir.tar.bz2"
}

if [ ! -e "$1" ]; then
	usage && exit
fi

if [ -d $1 ]
then
    tar cf $1.tar.bz2 --use-compress-prog=$HOME/bin/pbzip2 $1
elif [ -f $1 ]
then
# In this case, it bypasses core autodetection and sets 16 cores
# (16 "virtual cores" by HyperThreading):
    $HOME/bin/pbzip2 -p16 $1
else
    usage
fi
