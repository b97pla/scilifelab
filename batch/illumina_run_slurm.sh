#!/bin/sh

ILLUMINA_RUN=`basename $2`

sed -e "s/RUN/$ILLUMINA_RUN/g" ~roman/bin/sbatch_job_tmpl.sh | sbatch
