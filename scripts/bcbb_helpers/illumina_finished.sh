#!/bin/sh
# Wrapper for BCBB's illumina cronjob, checks whether the remote storage host is up
# and launches the bcl->qseq->fastq conversions.

# Ensure it loads the appropiate virtualenv
source ~/.bashrc

# Detect whether the remote (storage|analysis) machine is up and if so, launch the analysis
STORE_HOST=`grep store_host: ~/config/post_process.yaml | cut -f2 -d" "`
ping -q -c1 -w3 $STORE_HOST >& /dev/null && python ~/opt/bcbb/nextgen/scripts/illumina_finished_msg.py ~/config/universe_wsgi.ini ~/config/transfer_info.yaml $*

