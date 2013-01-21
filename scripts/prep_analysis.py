import os
import subprocess
import glob
import sys
from bcbio.pipeline.config_loader import load_config
if len(sys.argv) < 4:
        print """
Will copy files to the production account
Will then prepare the analysis directories data, intermediate, etc.

Usage:

prep_analysis.py  <flowcell_id> <project_id> <config_file> <dry_run>

	<flowcell_id>	eg: 121207_AD1H19ACXX 
	<project_id> 	eg: J.Doe_00_01
	<config_file>	--
	<dry_run>	y/n
        """
        sys.exit()

flowcell_id     = sys.argv[1]
project_name	= sys.argv[2]
config_file     = sys.argv[3]

try:
        config   = load_config(config_file)
        analysis = config['project']['analysis']
	INBOX	 = config['project']['INBOX']
except:
	sys.exit('ERROR: problem loading paths fom config file')

project_account	= INBOX.split('/')[2]

if sys.argv[4].lower()=='y':
	subprocess.call(["casava_data_delivery.py",project_name, flowcell_id, project_account, '-d'])
elif sys.argv[4].lower()=='n':
	subprocess.call(["casava_data_delivery.py",project_name, flowcell_id, project_account])
else:
	sys.exit("Argument <dry_run> is badly formated. Should be 'y' or 'n'")


if not os.path.isdir(os.path.join(analysis,project_name)):
	os.makedirs(os.path.join(analysis,project_name))

if not os.path.isdir(os.path.join(analysis,project_name,'data')):
	os.makedirs(os.path.join(analysis,project_name,'data'))

if not os.path.isdir(os.path.join(analysis,project_name,'intermediate')):
        os.makedirs(os.path.join(analysis,project_name,'intermediate'))

if not os.path.isdir(os.path.join(analysis,project_name,'data',flowcell_id)):
        os.makedirs(os.path.join(analysis,project_name,'data',flowcell_id))

if not os.path.isdir(os.path.join(analysis,project_name,'intermediate',flowcell_id)):
        os.makedirs(os.path.join(analysis,project_name,'intermediate',flowcell_id))

flist = glob.glob(os.path.join(INBOX,project_name,'*',flowcell_id,'*'))

for file_path in flist:
	file = file_path.split('/')[-1]
	#os.symlink(file_path , os.path.join(analysis,project_name,'data',flowcell_id,file))
	os.rename(file_path , os.path.join(analysis,project_name,'data',flowcell_id,file))
