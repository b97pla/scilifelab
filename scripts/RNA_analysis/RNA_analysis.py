import os
import sys
import subprocess

if len(sys.argv) < 7:
        print """
Usage:
        stand in 'intermediate' and run

        RNA_analysis.py <project id> <bed file> <gtf file> <mail> <config_file> <run dir 1> <run dir 2> ... <run dir N>

Arguments:
        <run dir i>
                - The name of the directory with the tophat_out_* -dirs.
                This is typically the same as the run name, such as
                20120323A_hiseq2000, but can be any name. The name of
                the run dir will also be the name set as the 'run name' 
                in the report. You might want to change this in the 
                rst file if your run dir name doesn't have an appropriate
                'run name'.

                - If more than one directory is given, the script will merge 
                the bamfiles from the diferent directories and do the analysis 
                on the merged runs

        <project id>
                - eg: M.Muurinen_11_01a

        <bed file>
                - reference gene model in bed format. Used by Ever-Seq
                to get gene body coverage and read distribution.

        <gtf fie>
                - reference annotation in gtf format, used by cufflinks and HTseq

        <mail>
                - mail adress for SLURM messages

        <config_file>
                - post_process.yaml"""
	sys.exit()

command=[os.environ['HOME']+'/opt/scilifelab/scripts/RNA_analysis/RNA_analysis.sh', '-p', sys.argv[1], '-b', sys.argv[2], '-g', sys.argv[3], '-m', sys.argv[4], '-c', sys.argv[5]] + sys.argv[6:]
command=' '.join(command)
os.system(command)

