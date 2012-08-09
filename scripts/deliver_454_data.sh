#! /bin/bash
#
# File: deliver_454_data.sh
# Created: Wed Jun 20 15:07:49 2012
# $Id: $
#
# Copyright (C) 2012 by Per Unneberg
#
# Author: Per Unneberg
#
# Description:
# Deliver 454 data to customer
# Simply performs an rsync of an input directory
# Logs are stored in /proj/a2010001/private/delivery_logs


LOGOUT=/proj/a2010001/private/delivery_logs/log.out
LOGERR=/proj/a2010001/private/delivery_logs/log.err

function usage {
    echo "Usage: $0 [-n] <src_dir> <dest_dir>"
    echo "Example: deliver_454_data.sh 2000_00_00_trolle /proj/x2000000/INBOX/J.Doe_00_01 "
}

function log_begin {
    echo "----------------------------"
    echo Running rsync $rsyncopts $1 $2
    echo `date`
    echo "----------------------------"
}
function log_end {
    echo "----------------------------"
    echo Finished
    echo `date`
    echo "----------------------------"
}


rsyncopts=-av
dry_run=0

while getopts ":n" opt; do
  case $opt in
    n)
      rsyncopts=-anv
      dry_run=1
      shift $((OPTIND-1))
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage && exit
      ;;
  esac
done


if [ ! -e "$1" ]; then
    usage && exit
fi

# dry run
if [ $dry_run == 1 ]; then
    rsync $rsyncopts $1 $2
    exit
fi

# Log start
log_begin
log_begin >> $LOGOUT

# do it
rsync $rsyncopts $1 $2 >> $LOGOUT 2>> $LOGERR

# End
log_end
log_end >> $LOGOUT
