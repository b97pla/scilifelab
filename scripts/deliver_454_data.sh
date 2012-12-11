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
    echo "Usage: $0 [-n12] <src_dir> <dest_dir>"
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

while getopts ":n12" opt; do
    case $opt in
	n)
	    rsyncopts="${rsyncopts}n"
	    dry_run=1
	    shift $((OPTIND-1))
	    ;;
	1)
	    rsyncopts="${rsyncopts}m --include='*/' --include='*01.sff' --include='1.TCA*' --include='1.cwf'   --exclude='*'"
	    shift $((OPTIND-1))
	    ;;
	2)
	    rsyncopts="${rsyncopts}m --include='*/' --include='*02.sff' --include='2.TCA*' --include='2.cwf'   --exclude='*'"
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
    echo rsync $rsyncopts $1 $2
    eval "rsync $rsyncopts $1 $2"
    exit
fi

# Log start
log_begin
log_begin >> $LOGOUT

# do it
eval "rsync $rsyncopts $1 $2 >> $LOGOUT 2>> $LOGERR"

# End
log_end
log_end >> $LOGOUT
