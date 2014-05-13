#!/usr/bin/env python

"""A module with lims help functions.

Maya Brandi, Science for Life Laboratory, Stockholm, Sweden.
"""

from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD
lims = Lims(BASEURI, USERNAME, PASSWORD)

"""process category dictionaries

In the lims_utils context, processes are categorised into groups that define, 
or are used to define a certain type of statusdb key. The categories and their 
processes are defined here:""" 'hh'

INITALQCFINISHEDLIB = {'24' : 'Customer Gel QC',
    '62' : 'qPCR QC (Library Validation) 4.0',
    '64' : 'Quant-iT QC (Library Validation) 4.0',
    '67' : 'Qubit QC (Library Validation) 4.0',
    '20' : 'CaliperGX QC (DNA)',
    '17' : 'Bioanalyzer QC (Library Validation) 4.0'}
INITALQC ={'63' : 'Quant-iT QC (DNA) 4.0',
    '65' : 'Quant-iT QC (RNA) 4.0',
    '66' : 'Qubit QC (DNA) 4.0',
    '68' : 'Qubit QC (RNA) 4.0',
    '24' : 'Customer Gel QC',
    '20' : 'CaliperGX QC (DNA)',
    '16' : 'Bioanalyzer QC (DNA) 4.0',
    '18' : 'Bioanalyzer QC (RNA) 4.0',
    '116' : 'CaliperGX QC (RNA)',
    '504' : 'Volume Measurement QC'}
AGRINITQC = {'7' : 'Aggregate QC (DNA) 4.0',
    '9' : 'Aggregate QC (RNA) 4.0'}
PREPREPSTART = {'74': 'Shear DNA (SS XT) 4.0',
    '304' : "Ligate 3' adapters (TruSeq small RNA) 1.0"}
PREPSTART = {'10' : 'Aliquot Libraries for Hybridization (SS XT)',
    '47' : 'mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0',
    '33' : 'Fragment DNA (TruSeq DNA) 4.0',
    '407' : 'Fragment DNA (Thruplex)',
    '308': 'Library Pooling (TruSeq Small RNA) 1.0',
    '117' : 'Applications Generic Process'}
PREPEND = {'157': 'Applications Finish Prep',
    '109' : 'CA Purification',
    '456' : 'Purification (ThruPlex)',
    '111' : 'Amplify Captured Libraries to Add Index Tags (SS XT) 4.0',
    '406' : 'End repair, size selection, A-tailing and adapter ligation (TruSeq PCR-free DNA) 4.0',
    '311': 'Sample Placement (Size Selection)'}
LIBVAL = {'62' : 'qPCR QC (Library Validation) 4.0',
    '64' : 'Quant-iT QC (Library Validation) 4.0',
    '67' : 'Qubit QC (Library Validation) 4.0',
    '20' : 'CaliperGX QC (DNA)',
    '17' : 'Bioanalyzer QC (Library Validation) 4.0',
    '24' : 'Customer Gel QC'}
AGRLIBVAL = {'8': 'Aggregate QC (Library Validation) 4.0'}
SEQSTART = {'23':'Cluster Generation (Illumina SBS) 4.0',
    '26':'Denature, Dilute and Load Sample (MiSeq) 4.0'}
DILSTART = {'40' : 'Library Normalization (MiSeq) 4.0',
    '39' : 'Library Normalization (Illumina SBS) 4.0'}
SEQUENCING = {'38' : 'Illumina Sequencing (Illumina SBS) 4.0',
    '46' : 'MiSeq Run (MiSeq) 4.0'}
WORKSET = {'204' : 'Setup Workset/Plate'}
SUMMARY = {'356' : 'Project Summary 1.3'}

PROJ_UDF_EXCEPTIONS = ['customer_reference','uppnex_id','reference_genome','application']

SAMP_UDF_EXCEPTIONS = ['customer_name','reads_requested_(millions)','min_reads','m_reads','dup_rm','status_auto','status_manual','average_size_bp','incoming_qc_status']

def get_udfs(udf_key, obj, udfs, exeptions = []):
    """Transforms udf names to statusdb keys (underscore and lowercase) and places them under
    details in obj. If exeptions are pased as argument, those will be placed on the 
    top level of obj

    Arguments:
    udf_key     string. name of key under wich udfs are collected.
    obj         dictionary. Eg instance of the Samples or Project classes
    udfs        udf dictionary
    exeptions   list of exception udf keys (underscore and lowercase)"""
    if not obj.has_key(udf_key):
        obj[udf_key]={}
    for key,val in udfs:
        try:
            val=_to_unicode(_from_unicode(val))
        except:
            pass
        db_key = key.replace(' ','_').lower()
        try:
            val = val.isoformat()
        except:
            pass
        if db_key in exeptions:
            obj[db_key] = val
        else:
            obj[udf_key][db_key] = val
    return obj



def get_sequencing_info(fc):
    """Input: a process object 'fc', of type 'Illumina Sequencing (Illumina SBS) 4.0',
    Output: A dictionary where keys are lanes 1,2,...,8, and values are lane artifact udfs"""
    fc_summary={}
    for iom in fc.input_output_maps:
        art = Artifact(lims,id = iom[0]['limsid'])
        lane = art.location[1].split(':')[0]
        if not fc_summary.has_key(lane):
            fc_summary[lane]= dict(art.udf.items()) #"%.2f" % val ----round??
            fc_summary[lane]['qc'] = art.qc_flag
    return fc_summary


def make_sample_artifact_maps(sample_name):
    """outin: connects each out_art for a specific sample to its 
    corresponding in_art and process. one-one relation"""
    
    outin = {}
    artifacts = lims.get_artifacts(sample_name = sample_name, type = 'Analyte') 
    for outart in artifacts:
        try:
            pro = outart.parent_process
            inarts = outart.input_artifact_list()
            for inart in inarts:
                for samp in inart.samples:
                    if samp.name == sample_name:
                        outin[outart.id] = (pro, inart.id)
        except:
            pass
    return outin

def get_analyte_hist_sorted(out_analyte, outin, inart = None):
    """Makes a history map of an analyte, using the outin-map 
    of the corresponding sample.
    The outin object is built up from analytes. This means that it will not 
    contain output-inpit info for processes wich have only files as output. 
    This is cusial since the outinobject is used for building upp the ANALYTE 
    history of a sample. If you want to make the analyte history based on a 
    resultfile, that is; if you want to give a resultfile as out_analyte here, 
    and be given the historylist of analytes and processes for that file, you 
    will also have to give the input artifact for the process that generated 
    the resultfile for wich you want to get the history. In other words, if you 
    want to get the History of the folowing scenario:        

    History --- > Input_analyte -> Process -> Output_result_file
    
    then the arguments to this function should be:
    out_analyte = Output_result_file
    inart = Input_analyte

    If you instead want the History of the folowing scenario:
    
    History --- > Input_analyte -> Process -> Output_analyte

    the you can skip the inart argument and only set:
    out_analyte = Output_analyte 
    """
    history = {}
    hist_list = []
    if inart:
        history, out_analyte = add_out_art_process_conection_list(inart, 
                                                    out_analyte, history)
        hist_list.append(inart)
    while outin.has_key(out_analyte):
        inart = outin[out_analyte][1]
        hist_list.append(inart)
        history, out_analyte = add_out_art_process_conection_list(inart, 
                                                        out_analyte, history)
    return history, hist_list

def add_out_art_process_conection_list(inart, out_analyte, history = {}):
    """This function populates the history dict with process info per artifact.
    Maps an artifact to all the processes where its used as input and adds this 
    info to the history dict. Obseve that the output artifavt for the input 
    atrifact in the historychain is given as input to this funktion. All 
    processes that the input artifakt has been involved in, but that are not 
    part of the historychain get the outart set to None. This is verry important."""
    processes = lims.get_processes(inputartifactlimsid = inart)
    for process in processes:
        outputs = map(lambda a: a.id, process.all_outputs())
        outart = out_analyte if out_analyte in outputs else None 
        step_info = {'date' : process.date_run,
                     'id' : process.id,
                     'outart' : outart,
                     'inart' : inart,
                     'type' : process.type.id,
                     'name' : process.type.name}
        if history.has_key(inart):
            history[inart][process.id] = step_info
        else:
            history[inart] = {process.id : step_info}
    return history, inart
