#!/usr/bin/env python
from uuid import uuid4
import time
from  datetime  import  datetime
import couchdb
#Make it backwards compatible
try:
    import bcbio.pipeline.config_utils as cl
except ImportError:
    import bcbio.pipeline.config_loader as cl

def load_couch_server(config_file):
    """loads couch server with settings specified in 'config_file'"""
    try:
        db_conf = cl.load_config(config_file)['statusdb']
        url = db_conf['username']+':'+db_conf['password']+'@'+db_conf['url']+':'+str(db_conf['port'])
        couch = couchdb.Server("http://" + url)
        return couch
    except KeyError:
        raise RuntimeError("\"statusdb\" section missing from configuration file.")

def find_or_make_key(key):
    if not key:
        key = uuid4().hex
    return key

def save_couchdb_obj(db, obj):
    """Updates ocr creates the object obj in database db."""
    dbobj = db.get(obj['_id'])
    time_log = datetime.utcnow().isoformat() + "Z"
    if dbobj is None:
        obj["creation_time"] = time_log
        obj["modification_time"] = time_log
        db.save(obj)
        return 'created'
    else:
        obj["_rev"] = dbobj.get("_rev")
        obj["modification_time"] = time_log
        dbobj["modification_time"] = time_log
        obj["creation_time"] = dbobj["creation_time"]
        if not comp_obj(obj, dbobj):
            db.save(obj)
            return 'uppdated'
    return 'not uppdated'

def save_couchdb_ref_obj(db, obj):
    """Updates ocr creates the object obj in database db."""
    dbobj = db.get(obj['_id'])
    time_log = datetime.utcnow().isoformat() + "Z"
    if dbobj is None:
        db.save(obj)
        return 'created'
    else:
        obj["_rev"] = dbobj.get("_rev")
        if not comp_obj(obj, dbobj):
            db.save(obj)
            return 'uppdated'
    return 'not uppdated'

def comp_obj(obj, dbobj):
    ####temporary
    if dbobj.has_key('entity_type'):
        if dbobj['entity_type']=='project_summary':
            obj=dont_load_status_if_20158_not_found(obj, dbobj)
    ###end temporary
    """compares the two dictionaries obj and dbobj"""
    keys = list(set(obj.keys() + dbobj.keys()))
    for key in keys:
        if (obj.has_key(key)) and dbobj.has_key(key):
            if (obj[key] != dbobj[key]):
                return False
        else:
            return False
    return True


def dont_load_status_if_20158_not_found(obj, dbobj):
    """compares the two dictionaries obj and dbobj"""
    if obj.has_key('samples') and dbobj.has_key('samples'):
        keys = list(set(obj['samples'].keys() + dbobj['samples'].keys()))
        for key in keys:
            if obj['samples'].has_key(key) and dbobj['samples'].has_key(key):
                if obj['samples'][key].has_key('status'):
                    if obj['samples'][key]['status'] == 'doc_not_found':
                        if dbobj['samples'][key].has_key('status'):
                            obj['samples'][key]['status'] = dbobj['samples'][key]['status']
                if obj['samples'][key].has_key('m_reads_sequenced'):
                    if obj['samples'][key]['m_reads_sequenced'] == 'doc_not_found':
                        if dbobj['samples'][key].has_key('m_reads_sequenced'):
                            obj['samples'][key]['m_reads_sequenced'] = dbobj['samples'][key]['m_reads_sequenced']
            try:
                if (obj['samples'][key]['status'] == 'doc_not_found') or (obj['samples'][key]['status'] == None):
                    obj['samples'][key].pop('status')
            except: pass
            try:
                if (obj['samples'][key]['m_reads_sequenced'] == 'doc_not_found') or (obj['samples'][key]['m_reads_sequenced'] == None):
                    obj['samples'][key].pop('m_reads_sequenced')
            except: pass
    return obj

def find_proj_from_view(proj_db, project_name):
    view = proj_db.view('project/project_name')
    for proj in view:
        if proj.key == project_name:
            return proj.value
    return None

def find_proj_from_samp(proj_db, sample_name):
    view = proj_db.view('samples/sample_project_name')
    for samp in view:
        if samp.key == sample_name:
            return samp.value
    return None

def find_samp_from_view(samp_db, proj_name):
    view = samp_db.view('names/id_to_proj')
    samps = {}
    for doc in view:
        if (doc.value[0] == proj_name) or (doc.value[0] == proj_name.lower()):
            samps[doc.key] = doc.value[1:3]
    return samps

def find_flowcell_from_view(flowcell_db, flowcell_name):
    view = flowcell_db.view('names/id_to_name')
    for doc in view:
        if doc.value:
            id = doc.value.split('_')[1]
            if (id == flowcell_name):
                return doc.key

def find_sample_run_id_from_view(samp_db,sample_run):
    view = samp_db.view('names/id_to_name')
    for doc in view:
        if doc.value == sample_run:
            return doc.key
    return None


