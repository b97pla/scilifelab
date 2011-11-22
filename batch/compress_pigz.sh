#!/bin/bash
#SBATCH -p node
#SBATCH -N 1
<<<<<<< HEAD
#SBATCH -t 5:00:00
#SBATCH -J pigz
#SBATCH -A a2010002
#SBATCH --mail-user=pontus.larsson@scilifelab.se
=======
#SBATCH -t 10:00:00
#SBATCH -J pbzip2
#SBATCH -A a2010002
#SBATCH --mail-user=roman.valls.guimera@scilifelab.se
>>>>>>> c9e32f980979c53142245d06367f6efc2b05d485
#SBATCH --mail-type=ALL

# Compress everything:
# http://www.dataspora.com/2010/08/the-seven-secrets-of-successful-data-scientists/
# ... using all cores with pbzip2

function usage {
    echo "Usage: $0 <-d> <-s> [dir|file] <dest_dir>"
    echo "Example: sbatch compress_slurm.sh bigfolder /store/longterm/mount/"
    echo "Output: file.bz2 if it's a single file or dir.tar.bz2, -d for decompressing"
}

<<<<<<< HEAD
if [ $1 == "-d" ]; then
	pigz -d $2
	exit
fi

=======
>>>>>>> c9e32f980979c53142245d06367f6efc2b05d485
if [ ! -e "$1" ]; then
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
<<<<<<< HEAD
    tar cf $1.tar.bz2 --use-compress-prog=pigz $1
=======
    tar cf $2`basename $1`.tar.bz2 --use-compress-prog=pbzip2 $1
>>>>>>> c9e32f980979c53142245d06367f6efc2b05d485
elif [ -f $1 ]
then
# In this case, it bypasses core autodetection and sets 16 cores
# (16 "virtual cores" by HyperThreading):
    pigz -p16 $1
else
    usage
fi
