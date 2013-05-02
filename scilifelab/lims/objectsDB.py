import codecs
from pprint import pprint
from statusDB_utils import *
from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD
lims = Lims(BASEURI, USERNAME, PASSWORD)
import os
import couchdb
import bcbio.pipeline.config_loader as cl
config_file = os.path.join(os.environ['HOME'], 'opt/config/post_process.yaml')
db_conf = cl.load_config(config_file)['couch_db']        
url = db_conf['maggie_login']+':'+db_conf['maggie_pass']+'@'+db_conf['maggie_url']+':'+str(db_conf['maggie_port']) 
samp_db = couchdb.Server("http://" + url)['samples']

class Lims2DB():
	def get_sample_status():
		"""ongoing,passed,aborted"""

	def get_preps(self, sample_name):
	        """Get preps and prep names; A,B,C... based on prep 
	        dates for sample_name. Returnes dict where keys are
	        prep_art_id and values are prep names."""
	        dates = {}
	        preps = {}
	        artifacts = lims.get_artifacts(sample_name = sample_name, process_type = 'mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0')
	        preps_keys = map(chr, range(65, 65+len(artifacts)))
	        for art in artifacts:
	                dates[art.id] = art.parent_process.date_run
	        for i, key in enumerate(sorted(dates,key= lambda x : dates[x])):
	                preps[key] = preps_keys[i]
	        return preps
	
	def get_prep_status(self, sample_name, prep_art_id):
	        """prep_art_id, artifact instance with parent process: 
	        mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0"""
	        status = {}
		average_size_bp = {}
	        artifacts = lims.get_artifacts(sample_name = sample_name, process_type = 'CA Purification')
	        for art in artifacts:
	                history, id_list = self.get_analyte_hist(art, sample_name)
	                if prep_art_id in id_list:
				average_size_bp[art.id] = dict(art.udf.items())['Size (bp)']
	                        status[art.id] = [art.qc_flag, art.parent_process.date_run]
	        return status, average_size_bp
	
	def get_analyte_hist(self, analyte, sample_name):
	        """List of dicts showing the history of an analyte.
	        sample_name has to be given since the analyte can
	        be a pool of many samples."""
	        history = []
	        id_list = []
	        while analyte:
	                step = {'out': analyte}
	                id_list.append(analyte.id)
	                try:
	                        step['process'] = analyte.parent_process
	                except:
	                        pass
	                inarts = analyte.input_artifact_list()
	                analyte = None
	                if len(inarts) > 0:
	                        for id in inarts:
	                                inart = Artifact(lims, id = id)
	                                if inart.samples[0].name == sample_name:
	                                        step['in'] = inart
	                                        analyte = inart
	                history.append(step)
	        return history, id_list

	def make_srm_id(self, lane_art, sample_name):
	        """Get sample_run_metrics id for sample: sample_name
	        run on lane: lane_art. Where lane_art is an artifact 
	        instance with parent process type: Cluster Generation 
	        (Illumina SBS) 4.0"""
	        try:
	                history, id_list = self.get_analyte_hist(lane_art, sample_name)
	                for step in history:
	                        if step.has_key('process'):
	                                if step['process'].type.name == 'Library Normalization (Illumina SBS) 4.0':
	                                        BARCODE = self.get_barcode(step['out'].reagent_labels[0])
	                                        DATE = step['process'].date_run[2:10].replace('-','')
	                FCID = lane_art.location[0].name
	                LANE = lane_art.location[1].split(':')[0]
	                return '_'.join([LANE,DATE,FCID,BARCODE])
	        except:
	                return None

	def get_barcode(self, name):
	        return name.split('(')[1].strip(')')

	def get_sample_runs(self, sample_name):
	        runs = {}
	        lane_artifacts = lims.get_artifacts(sample_name = sample_name, process_type = 'Cluster Generation (Illumina SBS) 4.0',type='Analyte')
	        for art in lane_artifacts:
	                runs[art.id] = self.make_srm_id(art, sample_name)
	        return runs


class ProjectDB(Lims2DB):
        """Convert project-udf-fields to project-couchdb-fields"""
        def __init__(self, project_id):
		self.lims_project = Project(lims,id = project_id)
		self.obj={}	# change name project_db_info
                self.udf_field_conv={'Name':'name',
                        'Queued':'queued',
                        'Portal ID':'Portal_id',
                        'Reference genome':'reference_genome',
                        'Sample type':'sample_type',
                        'Application':'application',
                        'Sequence units ordered (lanes)':'sequence_units_ordered_(lanes)',
                        'Sequencing platform':'sequencing_platform',
                        'Sequencing setup':'sequencing_setup',
                        'Library construction method':'library_construction_method',
                        'Bioinformatics':'bioinformatics',
                        'Disposal of any remaining samples':'disposal_of_any_remaining_samples',
                        'Type of project':'type_of_project',
                        'Invoice Reference':'invoice_reference',
                        'Uppmax Project':'uppnex_id',
                        'Uppmax Project Owner':'uppmax_project_owner',
                        'Custom Capture Design ID':'custom_capture_design_id',
                        'Customer project reference':'customer_reference',
                        'Customer Project Description':'customer_project_description',
                        'Project Comment':'project_comment',
                        'Delivery Report':'delivery_report'}
                for key, val in self.lims_project.udf.items():
                        self.obj[self.udf_field_conv[key]] = val
		samples = lims.get_samples(projectlimsid=self.lims_project.id)
		if len(samples)>0:
			self.obj['samples']={}
			for samp in samples:
				sampDB = SampleDB(samp.id) 
				self.obj['samples'][sampDB.name] = sampDB.obj


class SampleDB(Lims2DB):
        def __init__(self, sample_id):
                self.lims_sample = Sample(lims, id = sample_id)
                self.name = self.lims_sample.name
                self.obj={}
                self.udf_field_conv={'Name':'name',
                        'Customer Sample Name':'customer_name',
                        'Progress':'progress',
                        'Sequencing Method':'sequencing_method',
                        'Sequencing Coverage':'sequencing_coverage',
                        'Sample Type':'sample_type',
                        'Reference Genome':'reference_genome',
                        'Pooling':'pooling',
                        'Application':'application',
                        'Read Length':'read_length',
                        'Control?':'control?',
                        'Sample Buffer':'sample_buffer',
                        'Units':'units',
                        'Customer Volume':'customer_volume',
                        'Color':'color',
                        'Customer Conc.':'customer_conc',
                        'Customer Amount (ug)':'customer_amount_(ug)',
                        'Customer A260:280':'customer_A260:280',
                        'Conc Method':'conc_method',
                        'QC Method':'qc_method',
                        'Extraction Method':'extraction_method',
                        'Customer RIN':'customer_rin',
                        'Sample Links':'Sample Links',
                        'Sample Link Type':'Sample Link Type',
                        'Tumor Purity':'tumor_purity',
                        'Lanes Requested':'lanes_requested',
                        'Customer nM':'customer_nM',
                        'Customer Average Fragment Length':'customer_average_fragment_length',
                        'Reads Requested (millions)':'reads_requested_(millions)',
                        'Insert Size':'average_size_bp',
                        'Passed Initial QC':'incoming_QC_status',
                        'Passed Library QC':'prep_status',
                        '-DISCONTINUED-SciLifeLab ID':'sciLifeLab_ID',
			'-DISCONTINUED-Volume Remaining':'volume_remaining'}
                for key, val in self.lims_sample.udf.items():
                        self.obj[self.udf_field_conv[key]] = val
                runs = self.get_sample_runs(self.name)
                preps = self.get_preps(self.name)
                if len(preps)>0:
                        self.obj['library_prep'] = {}
                        for id, prep in preps.items():
                                status, size_bp  = self.get_prep_status(self.name, id)
                                self.obj['library_prep'][prep] = {'lims_id' : id,'prep_status' : status, 'average_size_bp': size_bp }
                                for run_id, sample_run_metrics in runs.items():
                                        history, id_list = self.get_analyte_hist(Artifact(lims,id = run_id), self.name)
                                        if id in id_list:
						sample_run_key = find_sample_run_id_from_view = (samp_db, sample_run_metrics)
                                                try:
                                                        self.obj['library_prep'][prep]['sample_run_metrics'][sample_run_metrics] = sample_run_key
                                                except:
                                                        self.obj['library_prep'][prep]['sample_run_metrics'] = {sample_run_metrics : sample_run_key }



