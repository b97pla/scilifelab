#! /bin/sh

#SBATCH -A a2010002
#SBATCH -t 6:00:00
#SBATCH -p core
#SBATCH -J md5sum_check
#SBATCH --qos=seqver

FCDIR=$1
LOG=`basename ${1}`.md5_check
echo $LOG
rm -f $LOG

cd $FCDIR
for d in `ls -d Unaligned*`
do
  find $d -type f -name "*.md5" -exec md5sum -c '{}' \; >> $LOG
  for f in `find $d -type f -name "*.fastq.gz"`
  do
    if [ ! -e ${f}.md5 ]
    then
      echo "FAIL: Missing .md5 file for $f" >> $LOG
    fi
  done
done
cd $OLDPWD
