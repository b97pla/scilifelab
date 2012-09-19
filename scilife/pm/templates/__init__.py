"""
Templates configuration
"""

import os
import sys
from mako.template import Template
from mako.lookup import TemplateLookup

TEMPLATE_ROOT = os.path.join(os.path.abspath(os.path.dirname(__file__)), "tpl")
MAKE_TEMPLATE_DIR = os.path.join(TEMPLATE_ROOT, "make")

def get_make_templates():
    tmpl = {'Makefile' : Template(filename=os.path.join(MAKE_TEMPLATE_DIR, 'Makefile'))}
    return tmpl



##############################
## Various templates
##############################
SBATCH_TEMPLATE = Template('''\
#!/bin/bash -l
TMPDIR=/scratch/$SLURM_JOB_ID

#SBATCH -A ${project_id}
#SBATCH -t ${time}
#SBATCH -o ${jobname}.stdout
#SBATCH -e ${jobname}.stderr
#SBATCH -J ${jobname}
#SBATCH -D ${workdir}
#SBATCH -p ${partition}
#SBATCH -n ${cores}
#SBATCH --mail-type=${mail_type}
#SBATCH --mail-user=${mail_user}
<%
if (constraint):
    constraint_str = "#SBATCH -C " + constraint
else:
    constraint_str = ""
%>
${constraint_str}
${header}
${command_str}
${footer}
''')
