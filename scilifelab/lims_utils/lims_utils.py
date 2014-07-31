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
POOLING = {'42': "Library Pooling (Illumina SBS) 4.0",
    '43': "Library Pooling (MiSeq) 4.0",
    '44': "Library Pooling (TruSeq Amplicon) 4.0",
    '45': "Library Pooling (TruSeq Exome) 4.0",
    '58': "Pooling For Multiplexed Sequencing (SS XT) 4.0",
    '255': "Library Pooling (Finished Libraries) 4.0",
    '308': "Library Pooling (TruSeq Small RNA) 1.0",
    '404': "Pre-Pooling (Illumina SBS) 4.0",
    '506': "Pre-Pooling (MiSeq) 4.0",
    '508': "Applications Pre-Pooling"}
PREPSTART = {'10' : 'Aliquot Libraries for Hybridization (SS XT)',
    '47' : 'mRNA Purification, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0',
    '33' : 'Fragment DNA (TruSeq DNA) 4.0',
    '407' : 'Fragment DNA (ThruPlex)',
    '308': 'Library Pooling (TruSeq Small RNA) 1.0',
    '117' : 'Applications Generic Process',
    '405' : 'RiboZero depletion, Fragmentation & cDNA synthesis (TruSeq RNA) 4.0'}
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
    '17' : 'Bioanalyzer QC (Library Validation) 4.0'}
LIBVALFINISHEDLIB = {'62' : 'qPCR QC (Library Validation) 4.0',
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
DEMULTIPLEX={'666' : 'Bcl Conversion & Demultiplexing (Illumina SBS) 4.0'}

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

def procHistory(proc, samplename):
    """Quick wat to get the ids of parent processes from the given process, 
    while staying in a sample scope"""
    hist=[]
    artifacts = lims.get_artifacts(sample_name = samplename, type = 'Analyte')
    not_done=True
    starting_art=proc.input_per_sample(samplename)[0].id
    while not_done:
        not_done=False 
        for o in artifacts:
            if o.id == starting_art:
                if o.parent_process is None:
                    #flow control : if there is no parent process, we can stop iterating, we're done.
                    not_done=False
                    break #breaks the for artifacts, we are done anyway.
                else:
                    not_done=True #keep the loop running
                hist.append(o.parent_process.id)
                for i in o.parent_process.all_inputs():
                    if i in artifacts:
                        # while increment
                        starting_art=i.id
                            
                        break #break the for allinputs, if we found the right one
                break # breaks the for artifacts if we matched the current one
    return hist 

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
