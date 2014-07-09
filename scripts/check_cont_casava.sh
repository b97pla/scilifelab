 #!/bin/bash
# Wraper for checking fastq_screen results for a specific project and flowcell
# Assumes casava structure

if [ $# -ne 2 ]; then
  echo "Usage:
        check_cont_casava.sh <flowcell id> <project id>

        Arguments:
        <flowcell id>
                - eg: 120127_BD0H2HACXX
        <project id>
                - eg: M.Muurinen_11_01a"
  exit
fi

fcID=$1
project_id=$2

dir=`pwd`

cd  /gulo/proj_nobackup/a2010002/illumina

grep Human ${project_id}/*/${fcID}/fastq_screen/*.txt | sed 's/:/\//g' | cut -f 2,6 -d '/'
echo ''
grep Mouse ${project_id}/*/${fcID}/fastq_screen/*.txt | sed 's/:/\//g' | cut -f 2,6 -d '/'
echo ''
grep Ecoli ${project_id}/*/${fcID}/fastq_screen/*.txt | sed 's/:/\//g' | cut -f 2,6 -d '/'
echo ''
grep Spruce ${project_id}/*/${fcID}/fastq_screen/*.txt | sed 's/:/\//g' | cut -f 2,6 -d '/'
echo ''
grep PhiX ${project_id}/*/${fcID}/fastq_screen/*.txt | sed 's/:/\//g' | cut -f 2,6 -d '/'


cd $dir


