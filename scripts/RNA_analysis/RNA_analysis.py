import os
import sys
import subprocess
import couchdb
import bcbio.pipeline.config_loader as cl
import optparse
from datetime import datetime

def find_proj_from_view(proj_db, proj_id):
    """Getting proj_id from the project database in statusdb via view"""
    view = proj_db.view('project/project_name')
    for proj in view:
        if proj.key == proj_id:
            return proj.value
    return None


def main(args,mail,conffile,analysis,stranded):
    project = args[0]
    runs = args[1:]
    conf = cl.load_config(conffile)
    port = conf['tools']['port']
    username = conf['tools']['login']
    password = conf['tools']['pass']
    URL = username+':'+password+'@'+conf['tools']['url']
    extra_arg = "#SBATCH " + conf['sbatch']['extra_arg']
    couch = couchdb.Server("http://" + URL + ':' +str(port))
    proj_db = couch['projects']
    key = find_proj_from_view(proj_db, project)
    info = proj_db[key]
    reference_genome = info['reference_genome']
    RNA_analysis_settings = conf['custom_algorithms']['RNA-seq analysis']
    refpath = RNA_analysis_settings[reference_genome]['genomepath']
    gtfpath = RNA_analysis_settings[reference_genome]['gtfpath']
    bedpath = RNA_analysis_settings[reference_genome]['bedpath']
    today = str(datetime.today().isoformat()).replace('-','_').split('.')[0].replace(':','_')
    command=[os.environ['HOME']+'/opt/scilifelab/scripts/RNA_analysis/RNA_analysis.sh', '-p', project, '-b', bedpath, '-g', gtfpath, '-m', mail, '-c', conffile, '-e', '"'+extra_arg+'"' ,'-a', str(analysis),'-s' , str(stranded),'-d',today] + runs
    command=' '.join(command)
    os.system(command)



if __name__ == '__main__':
    usage = """
Stand in "intermediate" and run:

RNA_analysis.py [options] <project id> <run dir 1> <run dir 2> ... <run dir N>

Arguments:
    <run dir i>
        - The name of the directory with the tophat_out_* -dirs.
        This is typically the same as the run name, such as
        8_130828_AD2736ACXX

        - If more than one directory is given, the script will merge 
        the bamfiles from the diferent directories and do the analysis 
        on the merged runs in a new directory called merged_<todays date>.

    <project id>
        - eg: M.Muurinen_11_01a"""
    parser = optparse.OptionParser(usage)

    parser.add_option('-c', '--config', action="store", dest="conffile", default=os.path.expanduser("~/opt/config/post_process.yaml"),
    help="Specify config file (post_process.yaml)")
    parser.add_option('-m', '--mail', action="store", dest="mail", default='',
    help="Specify a mailing address for SLURM mail notifications")
    parser.add_option('-s', '--stranded', action="store_true", dest="stranded", default="False",
    help="Run tophat with --librarytype fr-firststranded option for strand-specific RNAseq.")
    parser.add_option('-a', '--analysis_on_all', action="store_true", dest="analysis", default="False",
    help="Run ht-seq and cufflinks on all samples. Default will be to run only on merged samples.")
    (opts, args) = parser.parse_args()

    main(args,opts.mail,opts.conffile,opts.analysis,opts.stranded)
