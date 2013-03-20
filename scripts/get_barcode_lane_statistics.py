import sys
import os
import optparse
import glob
import yaml
import couchdb
import bcbio.google
import bcbio.pipeline.config_loader as cl


##
from statusdb import *
fc_con = FlowcellRunMetricsConnection(dbname="flowcells", username="u",password="p",url="maggie.scilifelab.se"§c_con.get_barcode_lane_statistics('M.Uhlen_12_08','P352_173B_index1','130220_AD1T9RACXX')

##

def get_barcode_lane_statistics(project_id, sample_id,lane, stat):
        for samp in stat:
                if (samp['Project'].replace('_','').replace('.','') == project_id) and (samp['Sample ID']==sample_id) and (samp["Lane"]==lane):
                        return samp['Mean Quality Score (PF)'], samp['% of >= Q30 Bases (PF)']
        return None

def find_fc_from_view(flowcell_db, flow_cell):
        view = flowcell_db.view('names/Barcode_lane_stat')
        for fc in view:
                if fc.key == flow_cell:
                        return fc.value
        return None





def main(project_id, flow_cell):
if 'u'=='u':
        CREDENTIALS_FILE = os.path.join(os.environ['HOME'], 'opt/config/gdocs_credentials')
        CONFIG_FILE = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')
        CONFIG = cl.load_config(CONFIG_FILE)
        URL = CONFIG['couch_db']['maggie_url'] + ':' + str(CONFIG['couch_db']['maggie_port'])
        couch = couchdb.Server("http://" + URL)
        flowcell_db = couch['flowcells']

        stat=find_fc_from_view(flowcell_db, '130220_AD1T9RACXX')



if __name__ == '__main__':
        main(sys.argv[1], sys.argv[2])
