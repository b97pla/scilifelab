#!/usr/bin/env python
import os
import sys
import urllib
import urllib2
from xml.etree import ElementTree

from optparse import OptionParser

BASE_URL = "http://www.ncrna.org/frnadb/api/"

def main(base_url=BASE_URL, organism="Homo sapiens", so_name=["tRNA","snRNA","snoRNA","rRNA","miRNA","piRNA"]):
    base_url=BASE_URL
    organism="Homo sapiens"
    so_name=["tRNA","snRNA","snoRNA","rRNA","miRNA","scRNA","RNase_MRP_RNA","RNase_P_RNA","SRP_RNA","telomerase_RNA","vault_RNA","Y_RNA"]
    entry_ids = {}
    for id in search(base_url,organism,so_name):
        entry_ids[id] = 1
    for id in entry_ids.keys():
        print id
    for entry_id in entry_ids.keys():
        for entry in fetch_entry(base_url, entry_id):
            print ">%s\n%s" % ("|".join([entry['id'],entry['ontology'],entry['organism'],entry['description']]),entry['sequence'])
    
def search(base_url, organism, so_name):
    so = " ".join(so_name) 
    query = urllib.quote("search/org=\"%s\" and sos any \"%s\"" % (organism,so))
    result = _read("%s%s" % (base_url,query))
    
    return _parse_frnadb_search_result(result)

def fetch_entry(base_url, entry):
    url = "%s%s" % (base_url,urllib.quote("entry/%s" % entry))
    return _parse_frnadb_entry_result(_read(url))

def _read(url):
    return urllib2.urlopen(url).read()

def _parse_frnadb_entry_result(result):
    entry_l = []
    result_e = ElementTree.fromstring(result)
    entry_list_e = result_e.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}entry_list")
    if entry_list_e is not None:
        for entry_e in entry_list_e.findall("{http://www.ncrna.org/frnadb/XMLSchema/entry}entry"):
            id = entry_e.get("id")
            sequence = entry_e.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}sequence").text
            description = entry_e.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}description").text or ""
            
            if not id or not sequence:
                continue
                
            ontology = entry_e.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}ontology")
            if ontology is not None:
                seq_ontology = ontology.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}sequence_ontology")
                ontology = seq_ontology.get("name")
            else:
                ontology = ""
                
            organism = entry_e.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}organism")
            if organism is not None:
                taxonomy = organism.find("{http://www.ncrna.org/frnadb/XMLSchema/entry}ncbi_taxonomy")
                organism = taxonomy.get("name")
            else:
                organism = ""
                
            entry_l.append({'id': id, 'sequence': sequence, 'description': description, 'ontology': ontology, 'organism': organism})
            
    return entry_l

def _parse_frnadb_search_result(result):
    result_l = []
    result_e = ElementTree.fromstring(result)
    entry_list_e = result_e.find("{http://www.ncrna.org/frnadb/XMLSchema/search}entry_list")
    if entry_list_e is not None:
        next_url = entry_list_e.get("next","")
        for entry_e in entry_list_e.findall("{http://www.ncrna.org/frnadb/XMLSchema/search}entry"):
            result_l.append(entry_e.get("id",None))
        if len(next_url) > 0:
            result_l.extend(_parse_frnadb_search_result(_read(next_url)))
    return result_l

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", "--organism", dest="organism", action="append", default=None)
    parser.add_option("-s", "--so_name", dest="so_name", action="append", default=None)
    parser.add_option("-u", "--base_url", dest="base_url", default=None)
    options, args = parser.parse_args()
    
    main(options.base_url, options.organism, options.so_name)

