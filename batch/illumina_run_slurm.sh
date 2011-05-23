#!/bin/sh
# Re-writes placeholders on SBATCH job template and launches job

function usage {
    echo "Usage: cd \$ANALYSIS_DIR/RUN_ANALYSIS_DIR && $0 <RUN_DIR>"
}

if [ ! -e "$1" ]; then
	usage && exit
fi

# Arguments compatibility with automated_initial_analysis.py from
# bcbb pipeline: http://github.com/brainstorm/bcbb
let $2=$1

RUN_ANALYSIS_DIR=`basename $1`

sed -e "s/RUN/$RUN_ANALYSIS_DIR/g" ~/bin/sbatch_job_tmpl.sh | sbatch
