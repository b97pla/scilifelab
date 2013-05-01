#!/usr/bin/env python
import time
from  datetime  import  datetime


def save_couchdb_obj(db, obj):
    dbobj = db.get(obj['_id'])
    time_log = datetime.utcnow().isoformat() + "Z"
    if dbobj is None:
        obj["creation_time"] = time_log
        obj["modification_time"] = time_log
        db.save(obj)
        return 'Created'
    else:
        obj["_rev"] = dbobj.get("_rev")
        del dbobj["modification_time"]
        obj["creation_time"] = dbobj["creation_time"]
        if not comp_obj(obj, dbobj):
            obj["modification_time"] = time_log
            db.save(obj)
            return 'Uppdated'
    return None

def comp_obj(obj, dbobj):
        for key in dbobj:
                if (obj.has_key(key)):
                        if (obj[key] != dbobj[key]):
                             return False
                else:
                        return False
        return True

def find_proj_from_view(proj_db, project_name):
        view = proj_db.view('project/project_name')
        for proj in view:
                if proj.key == project_name:
                        return proj.value
        return None

def find_samp_from_view(samp_db, proj_name):
        view = samp_db.view('names/id_to_proj')
        samps = {}
        for doc in view:
                if (doc.value[0] == proj_name)|(doc.value[0] == proj_name.lower()):
                        samps[doc.key] = doc.value[1:3]
        return samps

def find_sample_run_id_from_view(samp_db, sample_run):
        view = samp_db.view('names/id_to_name')
        for doc in view:
                if doc.value == sample_run:
        		return doc.key
	return None
