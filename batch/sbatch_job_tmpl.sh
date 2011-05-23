#!/bin/bash
# BCBB pipeline SLURM batch system launcher for automated_initial_analysis.py having arguments:
# Usage:
#    automated_initial_analysis.py <YAML config file> <flow cell dir>
#                                  [<YAML run information>]
#
# Uses a node (ideally with all its cores) to process the automated_initial_analysis

CONFIG=$1
RUN=$2

function usage {
    echo "Usage: $0 <YAML config file> [flow cell dir]"
}

if [[ "$VARIABLE" =~ .*_.*_.*_.* ]]
then
    echo -n "First argument appears to be a flow cell dir, using default config "
    RUN=$CONFIG
    CONFIG="$HOME/config/post_process.yaml"
    echo $CONFIG
fi

#if [ -d $RUN ]; then
#	echo "$RUN is not a flow cell dir"
#	usage && exit
#fi

# If no YAML config file provided ($1), use default value
CONFIG=${CONFIG:-"$HOME/config/post_process.yaml"}
ANALYSIS_DIR=`grep base_dir ~/config/post_process.yaml | cut -d" " -f4`
MAILTO="YOURMAIL@example.com"
# The pipeline does not rely on this variable.
# This is just to prevent system's /tmp filling up
# if some program *does* rely on it.
TMPDIR=/scratch/$SLURM_JOB_ID

if [ ! -e "$ANALYSIS_DIR" ]; then
	echo "$ANALYSIS_DIR not found or not defined in $CONFIG"
        usage && exit
fi

## SLURM directives, need adjustment, the above variables do not propagate here unfortunately
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 150:00:00
#SBATCH -J $RUN
#SBATCH -A a2010002
#SBATCH --mail-user=$MAILTO
#SBATCH --mail-type=ALL
#SBATCH -o ${RUN}.out
#SBATCH -e ${RUN}.err


cd $ANALYSIS_DIR/$RUN && automated_initial_analysis.py $CONFIG $ANALYSIS_DIR/store/$RUN $ANALYSIS_DIR/store/$RUN/run_info.yaml

# If the run preliminar analysis has finished successfully, generate an all lanes merged PDF summary and send it via email
if [ -e "run_summary.yaml" ]
then
	gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all_samples_summary.pdf *-summary.pdf

	echo "The files for run $RUN are available on $ANALYSIS_DIR/$RUN" | mailx -s "\[BCBB pipeline\] All samples PDF summary for $RUN" \
	-a all_samples_summary.pdf $MAILTO 
fi
