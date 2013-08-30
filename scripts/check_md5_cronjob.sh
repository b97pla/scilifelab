#!/bin/bash

ARCHIVE_DIR=/proj/a2010002/archive
CHECK_MD5=$HOME/bin/check_md5.sh

for fc_dir in `ls -d $ARCHIVE_DIR/1*XX`
do
    if [ ! -e $fc_dir/*md5_check ]
    then
        sbatch $CHECK_MD5 $fc_dir
    fi
done
