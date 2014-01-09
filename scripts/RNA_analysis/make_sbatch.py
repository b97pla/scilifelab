import sys
from bcbio.pipeline.config_loader import load_config

if len(sys.argv) < 7:
        print """
Usage:

make_sbatch.py <project> <partition> <time> <jobname> <mail> <config_file> [dependensy_list_arg]
	project			eg: a2012042
	partition		node/core
	time			eg: 13:00:00
        """
        sys.exit()
try:
        config    = load_config(sys.argv[6])
        extra_arg = "#SBATCH "+config['sbatch']['extra_arg']
except:
        extra_arg=""

try:
    dep="#SBATCH "+sys.argv[7]
except:
    dep=""

f=open(sys.argv[4]+".sh",'w')
print >>f, """#!/bin/bash -l
#SBATCH -A {0}
#SBATCH -p {1}
#SBATCH -t {2}
#SBATCH -e {3}.err
#SBATCH -o {3}.out
#SBATCH -J {3}
#SBATCH --mail-type=ALL
#SBATCH --mail-user={4}
{5}
{6}
""".format(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],extra_arg,dep)

