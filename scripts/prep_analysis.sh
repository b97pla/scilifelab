#! /bin/bash

if [ $# -ne 2 ]; then
echo "
Will make a fake delivery to /proj/a2012043/INBOX/ 
Will then prepare the analysis directories in /proj/a2012043/private/nobackup/projects/

Usage:
        prep_analysis.sh <project id> <flowcell id>

Arguments:
	<project id>	eg: J.Doe_00_01 
	<flowcell id> 	eg: 121207_AD1H19ACXX
"
exit
fi


project_name=$1
flowcell_id=$2
echo project_name
echo flowcell_id
casava_data_delivery.py $project_name $flowcell_id 'a2012043'

if [[ ! -d /proj/a2012043/private/nobackup/projects/${project_name} ]]; then
	mkdir /proj/a2012043/private/nobackup/projects/${project_name}
fi
if [[ ! -d /proj/a2012043/private/nobackup/projects/${project_name}/data ]]; then
	mkdir /proj/a2012043/private/nobackup/projects/${project_name}/data
fi
if [[ ! -d /proj/a2012043/private/nobackup/projects/${project_name}/intermediate ]]; then
	mkdir /proj/a2012043/private/nobackup/projects/${project_name}/intermediate
fi
if [[ ! -d /proj/a2012043/private/nobackup/projects/${project_name}/intermediate/${flowcell_id} ]]; then
	mkdir /proj/a2012043/private/nobackup/projects/${project_name}/intermediate/${flowcell_id}
fi
if [[ ! -d /proj/a2012043/private/nobackup/projects/${project_name}/data/${flowcell_id} ]]; then
	mkdir /proj/a2012043/private/nobackup/projects/${project_name}/data/${flowcell_id}
fi
cd /proj/a2012043/private/nobackup/projects/${project_name}/data/${flowcell_id}
ln -s /proj/a2012043/INBOX/${project_name}/*/${flowcell_id}/* .


