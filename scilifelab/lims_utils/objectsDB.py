#!/usr/bin/env python

"""A module for building up the project objects that build up the project database on 
statusdb with lims as the main source of information.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""
import load_status_from_google_docs  ### Temorary solution untill 20158 implemented in LIMS!!!
import codecs
from scilifelab.google import _to_unicode, _from_unicode
from pprint import pprint
from genologics.lims import *
from lims_utils import *
from scilifelab.db.statusDB_utils import *
from helpers import *
import os
import couchdb
import bcbio.pipeline.config_utils as cl
import time
from datetime import date


###  Functions ###

def udf_dict(element, dict = {}):
    for key, val in element.udf.items():
        key = key.replace(' ', '_').lower().replace('.','')
        dict[key] = val
    return dict

def get_last_first(process_list, last=True):
    if process_list:
        process = process_list[0]
        for pro in process_list:
            new_date = int(pro['date'].replace('-',''))
            old_date = int(process['date'].replace('-',''))
            if last and (new_date > old_date):
                process = pro
            elif not last and (new_date < old_date):
                process = pro
        return process
    else:
        return None

### Classes  ###

class ProjectDB():
    """Instances of this class holds a dictionary formatted for building up the project database on statusdb. 
    Source of information come from different lims artifacts and processes. A detailed documentation of the 
    source of all values is found in: 
    https://docs.google.com/a/scilifelab.se/document/d/1OHRsSI9btaBU4Hb1TiqJ5wwdRqUQ4BAyjJR-Nn5qGHg/edit#"""
    def __init__(self, lims_instance, project_id, samp_db):
        self.lims = lims_instance 
        self.samp_db = samp_db
        self.lims_project = Project(self.lims,id = project_id)
        self.preps = ProcessInfo(self.lims , self.lims.get_processes(projectname = self.lims_project.name, type = AGRLIBVAL.values()))
        runs = self.lims.get_processes(projectname = self.lims_project.name, type = SEQUENCING.values())
        self.runs = ProcessInfo(self.lims, runs)
        project_summary = self.lims.get_processes(projectname = self.lims_project.name, type = SUMMARY.values())
        self.project = {'source' : 'lims',
            'application' : None,
            'samples':{},
            'open_date' : self.lims_project.open_date,
            'close_date' : self.lims_project.close_date,
            'entity_type' : 'project_summary',
            'contact' : self.lims_project.researcher.email,
            'project_name' : self.lims_project.name,
            'project_id' : self.lims_project.id}
        self.project = get_udfs('details', self.project, self.lims_project.udf.items(), PROJ_UDF_EXCEPTIONS)
        if dict(self.lims_project.researcher.lab.udf.items()).has_key('Affiliation'):
            self.project['affiliation'] = dict(self.lims_project.researcher.lab.udf.items())['Affiliation']
        if len(project_summary) == 1:
            self.project = get_udfs('project_summary', self.project, project_summary[0].udf.items())
        elif len(project_summary) > 1:
            print 'Warning. project summary process run more than once'

        #################Temporary solution untill 20158 implemented in lims >>>>>>>>>>>>>>>>>>>>>>>
        ## can be reooved when all project opened before 2014-01-27 have been closed for more than 60 days
        ## then we also need to block old projects so that they are not overwriten in case of someone manualy 
        ## updating it with the -p flagg
        opened = self.lims_project.open_date
        if opened:
            googledocs_status = {}
            try:
                googledocs_status = load_status_from_google_docs.get(self.lims_project.name)
            except:
                print 'issues finding status from 20158'
                pass
        ## Finish Date = last seq date if proj closed. Will be removed and feched from lims.
        seq_finished = None
        if self.lims_project.close_date and len(runs) > 0:
            d = '2000-10-10'
            for run in runs:
                try:
                    new_date = dict(run.udf.items())['Finish Date'].isoformat()
                    if comp_dates(d,new_date):
                        d = new_date
                    seq_finished = d
                except:
                    pass
        self.project['sequencing_finished'] = seq_finished
        #Temporary solution untill 20158 implemented in lims <<<<<<<<<<<<<<<<<<<<<<<

        ## Getting sample info
        samples = self.lims.get_samples(projectlimsid = self.lims_project.id)
        self.project['no_of_samples'] = len(samples)
        if len(samples) > 0:
            self.project['first_initial_qc'] = '3000-10-10'
            for samp in samples: 
                sampDB = SampleDB(self.lims,
                                samp.id,
                                self.project['project_name'],
                                self.samp_db,
                                self.project['application'],
                                self.preps.info,
                                self.runs.info,
                                googledocs_status) #googledocs_status Temporary solution untill 20158 implemented in lims!!
                self.project['samples'][sampDB.name] = sampDB.obj
##### initial qc fixa
                try:
                    initial_qc_start_date = self.project['samples'][sampDB.name]['initial_qc']['start_date']
                    if comp_dates(initial_qc_start_date,self.project['first_initial_qc']):
                        self.project['first_initial_qc'] = initial_qc_start_date
                except:
                    pass
        self.project = delete_Nones(self.project)

class ProcessInfo():
    """This class takes a list of process type names. Eg 
    'Aggregate QC (Library Validation) 4.0' and forms  a dict with info about 
    all processes of the type specified in runs which the project has gon through.

    info = {24-8460:{'finish_date':'2013-04-20', 
              'start_date',
              'run_id':'24-8460',
              'samples':{'P424_111':{in_art_id1 : [in_art1, out_art1],
                         in_art_id2: [in_art2, out_art2]},
                     'P424_115': ...},
                       ...},
        '24-8480':...}"""
    def __init__(self, lims_instance, runs):
        self.lims = lims_instance
        self.info = self.get_run_info(runs)
    def get_run_info(self, runs):
        run_info = {}
        for run in runs:
            run_info[run.id] = {'type' : run.type.name ,
                                'start_date': run.date_run,
                                'samples' : {}}
            run_udfs = dict(run.udf.items())
            try:
                run_info[run.id]['run_id'] = run_udfs["Run ID"]
            except:
                pass
            try:
                run_info[run.id]['finish_date'] = run_udfs['Finish Date'].isoformat()
            except:
                run_info[run.id]['finish_date'] = None
                pass
            in_arts=[]
            for IOM in run.input_output_maps:
                in_art_id = IOM[0]['limsid']
                in_art = Artifact(self.lims, id= in_art_id)
                out_art_id = IOM[1]['limsid']
                out_art = Artifact(self.lims, id= out_art_id)
                samples = in_art.samples
                if in_art_id not in in_arts:
                    in_arts.append(in_art_id)
                    for samp in samples:
                        if not samp.name in run_info[run.id]['samples'].keys():
                            run_info[run.id]['samples'][samp.name] = {}
                        run_info[run.id]['samples'][samp.name][in_art_id] = [in_art, out_art]
        return run_info


class SampleDB():
    """Instances of this class holds a dictionary formatted for building up the 
    samples in the project database on status db. Source of information come 
    from different lims artifacts and processes. A detailed documentation of the
     source of all values is found in
    https://docs.google.com/a/scilifelab.se/document/d/1OHRsSI9btaBU4Hb1TiqJ5wwdRqUQ4BAyjJR-Nn5qGHg/edit#"""
    def __init__(self,lims_instance , sample_id, project_name, samp_db, 
                        application = None, prep_info = [], run_info = [],
                        googledocs_status = {}): 
      # googledocs_status temporary solution untill 20158 implemented in lims!!
        self.lims = lims_instance
        self.samp_db = samp_db
        self.AgrLibQCs = prep_info
        self.lims_sample = Sample(self.lims, id = sample_id)
        self.name = self.lims_sample.name
        self.application = application
        self.outin = make_sample_artifact_maps(self.name)
        self.obj = get_udfs('details', {}, 
                                self.lims_sample.udf.items(), 
                                SAMP_UDF_EXCEPTIONS)
        self.obj['scilife_name'] = self.name
        self.obj['well_location'] = self.lims_sample.artifact.location[1]
        preps = self._get_preps_and_libval()
        if preps:
            runs = self.get_sample_run_metrics(run_info, preps)
            for prep_id in runs.keys():
                if preps.has_key(prep_id):
                    preps[prep_id]['sample_run_metrics'] = runs[prep_id]
            self.obj['library_prep'] = self._get_prep_leter(preps)
        initialqc = self._get_initialqc()
        self.obj['initial_qc'] = initialqc if initialqc else None
        if self.application in ['Finished library', 'Amplicon']:
            init_qc = INITALQCFINISHEDLIB.values()
        else:
            init_qc = INITALQC.values() 
        self.obj['first_initial_qc_start_date'] = self._get_firts_day(self.name,
                                                                        init_qc)
        self.obj['first_prep_start_date'] = self._get_firts_day(self.name, 
                                    PREPSTART.values() + PREPREPSTART.values())
        if googledocs_status and self.name in googledocs_status.keys():
            self.obj['status'] = googledocs_status[self.name][0]
            self.obj['m_reads_sequenced'] = googledocs_status[self.name][1]
        self.obj = delete_Nones(self.obj)

    def _get_firts_day(self, sample_name ,process_list, last_day = False):
        """process_list is a list of process type names, sample_name is a 
        sample name :)"""
        arts = self.lims.get_artifacts(sample_name = sample_name, 
                                        process_type = process_list)
        days = map(lambda a: a.parent_process.date_run , arts)
        days = filter(lambda d: d!=None  , days)
        if days:
            return max(days) if last_day else min(days)
        else:
            return None

    def get_barcode(self, reagent_label):
        """Extracts barcode from list of artifact.reagent_labels"""
        if reagent_label:
            try:
                index = reagent_label.split('(')[1].strip(')')
            except:
                index = reagent_label
        else:
            return None
        return index

    def get_sample_run_metrics(self, SeqRun_info, preps):
        """Input: SeqRun_info - instance of the ProcessInfo class with 
        SEQUENCING processes as argument
        For each SEQUENCING process run on the sample, this function steps 
        bacward in the artifact history of the input artifact of the SEQUENCING 
        process to find the folowing information:

        dillution_and_pooling_start_date  date-run of DILSTART step
        sequencing_start_date             date-run of SEQSTART step
        sequencing_run_QC_finished        date-run of SEQUENCING step
        sequencing_finish_date            udf ('Finish Date') of SEQUENCING step
        sample_run_metrics_id             The sample database (statusdb) _id for
                                          the sample_run_metrics corresponding 
                                           to the run, sample, lane in question.
        samp_run_met_id = lane_date_fcid_barcode            
            date and fcid:  from udf ('Run ID') of the SEQUENCING step. 
            barcode:        The reagent-lables of the input artifact of process 
                            type AGRLIBVAL
            lane:           from the location of the input artifact to the 
                            SEQUENCING step    
        preps are defined as the id of the PREPSTART step in the artifact 
        history. If appllication== Finished library, prep is defined as 
        "Finnished". These keys are used to connect the seqeuncing steps to the 
        correct preps."""
        sample_runs = {}
        for id, run in SeqRun_info.items():
            if run['samples'].has_key(self.name) and run.has_key('run_id'):
                date = run['run_id'].split('_')[0]
                fcid = run['run_id'].split('_')[3]
                run_type = run['type']
                for id , arts in run['samples'][self.name].items():
                    lane_art = arts[0]
                    outart = arts[1]
                    if run_type == "MiSeq Run (MiSeq) 4.0":
                        lane = lane_art.location[1].split(':')[1]
                    else:
                        lane = lane_art.location[1].split(':')[0]
                    hist_sort, hist_list = get_analyte_hist_sorted(outart.id, 
                                                        self.outin, lane_art.id)
                    steps = ProcessSpec(hist_sort, hist_list, self.application)
                    if self.application in ['Finished library', 'Amplicon']:
                        key = 'Finished'
                    elif steps.preprepstart:
                        key = steps.preprepstart['id']
                    elif steps.prepstart:
                        key = steps.prepstart['id'] 
                    else:
                        key = None 
                    if key:
                        if preps[key].has_key('reagent_label'):
                            barcode = self.get_barcode(preps[key]['reagent_label'])
                            samp_run_met_id = '_'.join([lane, date, fcid, barcode])
                            dict = {'sample_run_metrics_id':find_sample_run_id_from_view(self.samp_db, samp_run_met_id),
                                'dillution_and_pooling_start_date' : steps.dilstart['date'] if steps.dilstart else None,
                                'sequencing_start_date' : steps.seqstart['date'] if steps.seqstart else None,
                                'sequencing_run_QC_finished' : run['start_date'],
                                'sequencing_finish_date' : run['finish_date']}
                            dict = delete_Nones(dict)
                            if not sample_runs.has_key(key):
                                sample_runs[key] = {}
                            sample_runs[key][samp_run_met_id] = dict
        return sample_runs

    def _get_prep_leter(self, prep_info):
        """Get preps and prep names; A,B,C... based on prep dates for 
        sample_name. 
        Output: A dict where keys are prep_art_id and values are prep names."""
        dates = {}
        prep_info_new = {}
        preps_keys = map(chr, range(65, 65+len(prep_info)))
        if len(prep_info) == 1:
            prep_info_new['A'] = prep_info.values()[0]
        else:
            for key, val in prep_info.items():
                if val['pre_prep_start_date']:
                    dates[key] = val['pre_prep_start_date']
                else:
                    dates[key] = val['prep_start_date']
            for i, key in enumerate(sorted(dates,key= lambda x : dates[x])):
                prep_info_new[preps_keys[i]] = delete_Nones(prep_info[key])
        return prep_info_new

    def _get_preps_and_libval(self):
        """"""
        top_level_agrlibval_steps = self._get_top_level_agrlibval_steps()
        preps = {}
        very_last_libval_key = {}
        for AgrLibQC_id in top_level_agrlibval_steps.keys():
            AgrLibQC_info = self.AgrLibQCs[AgrLibQC_id]
            if AgrLibQC_info['samples'].has_key(self.name):
                inart, outart = AgrLibQC_info['samples'][self.name].items()[0][1]
                hist_sort, hist_list = get_analyte_hist_sorted(outart.id,
                                                               self.outin,
                                                               inart.id)
                steps = ProcessSpec(hist_sort, hist_list, self.application)
                prep = Prep()
                prep.set_prep_info(steps, self.application)
                if not preps.has_key(prep.id2AB) and prep.id2AB:
                    preps[prep.id2AB] = prep.prep_info
                if prep.pre_prep_library_validations and prep.id2AB:
                    preps[prep.id2AB]['pre_prep_library_validation'].update(prep.pre_prep_library_validations)
                if prep.library_validations and prep.id2AB:
                    preps[prep.id2AB]['library_validation'].update(prep.library_validations)
                    last_libval_key = max(prep.library_validations.keys())
                    last_libval = prep.library_validations[last_libval_key]
                    if not very_last_libval_key.has_key(prep.id2AB) or (last_libval_key > very_last_libval_key[prep.id2AB]):
                        very_last_libval_key[prep.id2AB] = last_libval_key
                        if last_libval.has_key('prep_status'):
                            preps[prep.id2AB]['prep_status'] = last_libval['prep_status']
                        preps[prep.id2AB]['reagent_label'] = self._pars_reagent_labels(steps, last_libval)
        if preps.has_key('Finished'):
            preps['Finished']['reagent_label'] = self.lims_sample.artifact.reagent_labels[0]
            preps['Finished'] = delete_Nones(preps['Finished'])
        
        return preps


    def _pars_reagent_labels(self, steps, last_libval):
        if steps.firstpoolstep:
            inart = Artifact(lims, id = steps.firstpoolstep['inart'])
            if len(inart.reagent_labels) == 1:
                return inart.reagent_labels[0]
        if last_libval.has_key('reagent_labels'): 
            if len(last_libval['reagent_labels']) == 1:
                return last_libval['reagent_labels'][0]
            return None
        return None

    def _get_initialqc(self):
        agr_qc = AGRINITQC
        outarts = self.lims.get_artifacts(sample_name = self.name, 
                                                process_type = agr_qc.values())
        parent_proc = map(lambda a: a.parent_process ,outarts)
        initialqc = {}
        if outarts:
            outart = Artifact(lims, id = max(map(lambda a: a.id, outarts)))
            latestInitQc = outart.parent_process
            inart = latestInitQc.input_per_sample(self.name)[0].id
            hist_sort, hist_list = get_analyte_hist_sorted(outart.id, 
                                                        self.outin, inart)
            if hist_list:
                iqc = InitialQC(hist_sort, hist_list)
                initialqc = delete_Nones(iqc.set_initialqc_info())
        return delete_Nones(initialqc)       

    def _get_top_level_agrlibval_steps(self):
        topLevel_AgrLibQC={}
        for AgrLibQC_id, AgrLibQC_info in self.AgrLibQCs.items():
            if AgrLibQC_info['samples'].has_key(self.name):
                topLevel_AgrLibQC[AgrLibQC_id]=[]
                inart, outart = AgrLibQC_info['samples'][self.name].items()[0][1]
                hist_sort, hist_list = get_analyte_hist_sorted(outart.id,
                                                          self.outin, inart.id)
                for inart in hist_list:
                    proc_info = hist_sort[inart]
                    proc_info = filter(lambda p : 
                             (p['type'] in AGRLIBVAL.keys()),proc_info.values())
                    
                    proc_ids = map(lambda p : p['id'], proc_info) 
                    topLevel_AgrLibQC[AgrLibQC_id] = topLevel_AgrLibQC[AgrLibQC_id] + proc_ids
        for AgrLibQC, LibQC in topLevel_AgrLibQC.items():
            LibQC=set(LibQC)
            if LibQC:
                for AgrLibQC_comp, LibQC_comp in topLevel_AgrLibQC.items():
                    if AgrLibQC_comp != AgrLibQC:
                        LibQC_comp=set(LibQC_comp)
                        if LibQC.issubset(LibQC_comp) and topLevel_AgrLibQC.has_key(AgrLibQC):
                            topLevel_AgrLibQC.pop(AgrLibQC)
        return topLevel_AgrLibQC

class InitialQC():
    """"""
    def __init__(self, hist_sort, hist_list, finnished_lib = False):
        self.init_qc = INITALQCFINISHEDLIB if finnished_lib else INITALQC
        self.agr_qc = AGRLIBVAL if finnished_lib else AGRINITQC
        self.initialqcend = None
        self.initialqcends = []
        self.initialqcs = []
        self.initialqstart = None
        self._set_initialqc_processes(hist_sort, hist_list)

    def _set_initialqc_processes(self, hist_sort, hist_list):
        for inart in hist_list:
            art_steps = hist_sort[inart]
            # INITALQCEND - get last agr initialqc val step after prepreplibval
            self.initialqcends += filter(lambda pro: pro['type'] in self.agr_qc, 
                                                            art_steps.values())
            # INITALQCSTART - get all lib val step after prepreplibval
            self.initialqcs += filter(lambda pro: pro['type'] in self.init_qc,
                                                            art_steps.values())
        self.initialqcend = get_last_first(self.initialqcends, last = True)
        self.initialqstart =  get_last_first(self.initialqcs, last = False)

    def set_initialqc_info(self):
        initialqc_info = {}
        initialqc_info['start_date'] = self.initialqstart['date'] if self.initialqstart else None
        if self.initialqcend:
            inart = Artifact(lims, id = self.initialqcend['inart'])
            process = Process(lims,id = self.initialqcend['id'])
            initialqc_info = udf_dict(inart, initialqc_info)
            initials = process.technician.initials
            initialqc_info['initials'] = initials
            initialqc_info['finish_date'] = self.initialqcend['date']
            initialqc_info['initial_qc_status'] = inart.qc_flag
        return initialqc_info


class ProcessSpec():
    def __init__(self, hist_sort, hist_list, application):
        self.application = application
        self.libvalends = []                
        self.libvalend = None               
        self.libvals = []
        self.libvalstart = None
        self.prepend = None                 
        self.prepstarts = []
        self.prepstart = None
        self.prepreplibvalends = []         
        self.prepreplibvalend = None        
        self.prepreplibvals = []
        self.prepreplibvalstart = None
        self.preprepstarts = []             
        self.prepends = []
        self.preprepstart = None
        self.workset = None                 
        self.worksets = []
        self.seqstarts = []
        self.seqstart = None
        self.dilstart = None
        self.dilstarts = []
        self.poolingsteps = []
        self.firstpoolstep = None
        self._set_prep_processes(hist_sort, hist_list)

    def _set_prep_processes(self, hist_sort, hist_list):
        hist_list.reverse()
        for inart in hist_list:
            prepreplibvalends = []
            art_steps = hist_sort[inart]
            # 1) PREPREPSTART
            self.preprepstarts += filter(lambda pro: (pro['type'] in 
                            PREPREPSTART and pro['outart']), art_steps.values()) #and pro['outart'] ##26 may

            if self.preprepstarts and not self.prepends: 
                # 2)PREPREPLIBVALSTART PREPREPLIBVALEND
                self.prepreplibvals += filter(lambda pro: (pro['type'] in 
                                                LIBVAL), art_steps.values())
                self.prepreplibvalends += filter(lambda pro: pro['type'] in
                                                AGRLIBVAL, art_steps.values())
            elif self.application in ['Finished library', 'Amplicon']: 
                # 6) LIBVALSTART LIBVALEND
                self.libvals += filter(lambda pro: pro['type'] in
                                                LIBVALFINISHEDLIB, art_steps.values())
                self.libvalends += filter(lambda pro: pro['type'] in
                                                AGRLIBVAL, art_steps.values())
            elif self.prepends: 
                # 6) LIBVALSTART LIBVALEND
                self.libvals += filter(lambda pro: pro['type'] in
                                                LIBVAL, art_steps.values())
                self.libvalends += filter(lambda pro: pro['type'] in
                                                AGRLIBVAL, art_steps.values())

            # 4) PREPSTART
            self.prepstarts += filter(lambda pro: (pro['type'] in 
                            PREPSTART) and pro['outart'], art_steps.values()) 
            # 5) PREPEND            - get latest prep end
            self.prepends += filter(lambda pro: (pro['type'] in 
                            PREPEND) and pro['outart'] , art_steps.values())

            # 8) WORKSET            - get latest workset
            self.worksets += filter(lambda pro: (pro['type'] in 
                            WORKSET) and pro['outart'], art_steps.values()) 

            # 9) SEQSTART dubbelkolla
            if not self.seqstarts:
                self.seqstarts = filter(lambda pro: (pro['type'] in SEQSTART) 
                                        and pro['outart'], art_steps.values())
            # 10) DILSTART dubbelkolla
            if not self.dilstarts:
                self.dilstarts = filter(lambda pro: (pro['type'] in DILSTART) 
                                        and pro['outart'], art_steps.values())
            # 11) POOLING STEPS
            self.poolingsteps += filter(lambda pro: (pro['type'] in
                                        POOLING), art_steps.values()) 
        self.workset = get_last_first(self.worksets) 
        self.libvalstart = get_last_first(self.libvals, last = False)
        self.libvalend = get_last_first(self.libvalends)
        self.prepreplibvalend = get_last_first(self.prepreplibvalends)
        self.prepstart = get_last_first(self.prepstarts, last = False)
        self.prepend = get_last_first(self.prepends)
        self.prepreplibvalstart = get_last_first(self.prepreplibvals, 
                                                            last = False)
        self.preprepstart = get_last_first(self.preprepstarts, last = False)
        self.firstpoolstep = get_last_first(self.poolingsteps, last = False)
        self.dilstart = get_last_first(self.dilstarts, last = False)
        self.seqstart = get_last_first(self.seqstarts, last = False)

class Prep():
    def __init__(self):
        self.prep_info = {
            'reagent_label': None,
            'library_validation':{},
            'pre_prep_library_validation':{},
            'prep_start_date': None,
            'prep_finished_date': None,
            'prep_id': None,
            'workset_setup': None,
            'pre_prep_start_date' : None}
        self.id2AB = None
        self.library_validations = {}
        self.pre_prep_library_validations = {}
        self.lib_val_templ = {
            'start_date' : None,
            'finish_date' : None,
            'well_location' : None,
            'prep_status' : None,
            'reagent_labels' : None,
            'average_size_bp' : None,
            'initials' : None}

    def set_prep_info(self, steps, aplication):
        if aplication in ['Amplicon', 'Finished library']:
            self.id2AB = 'Finished'
        else:
            if steps.prepstart:
                self.prep_info['prep_start_date'] = steps.prepstart['date']
            if steps.prepend:
                self.prep_info['prep_finished_date'] = steps.prepend['date']
                self.prep_info['prep_id'] = steps.prepend['id']
            if steps.workset:
                self.prep_info['workset_setup'] = steps.workset['id']
            if steps.preprepstart:
                self.prep_info['pre_prep_start_date'] = steps.preprepstart['date']
                self.id2AB = steps.preprepstart['id']
                if steps.preprepstart['outart']:
                    self.prep_info = udf_dict(Artifact(lims, 
                            id = steps.preprepstart['outart']), self.prep_info)
            elif steps.prepstart:
                self.id2AB = steps.prepstart['id']
                if steps.prepstart['outart']:
                    self.prep_info = udf_dict(Artifact(lims, 
                                id = steps.prepstart['outart']), self.prep_info)
        if steps.libvalend:
            self.library_validations = self._get_lib_val_info(steps.libvalends,
                                                            steps.libvalstart)
        if steps.prepreplibvalend:
            self.pre_prep_library_validations = self._get_lib_val_info(
                              steps.prepreplibvalends, steps.prepreplibvalstart)

        
    def _get_lib_val_info(self, agrlibQCsteps, libvalstart):
        library_validations = {}
        start_date = libvalstart['date'] if (libvalstart and 
                                         libvalstart.has_key('date')) else None
        for agrlibQCstep in agrlibQCsteps:
            library_validation = self.lib_val_templ
            inart = Artifact(lims, id = agrlibQCstep['inart'])
            if agrlibQCstep.has_key('date'):
                library_validation['finish_date'] = agrlibQCstep['date']
            library_validation['start_date'] = start_date
            library_validation['well_location'] = inart.location[1]
            library_validation['prep_status'] = inart.qc_flag
            library_validation['reagent_labels'] = inart.reagent_labels
            library_validation = udf_dict(inart, library_validation)
            initials = Process(lims, id = agrlibQCstep['id']).technician.initials
            if initials:
                library_validation['initials'] = initials
            if library_validation.has_key("size_(bp)"):
                average_size_bp = library_validation.pop("size_(bp)")
                library_validation["average_size_bp"] = average_size_bp
            library_validations[agrlibQCstep['id']] = delete_Nones(library_validation)
        return delete_Nones(library_validations) 
