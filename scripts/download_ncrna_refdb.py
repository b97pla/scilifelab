#!/usr/bin/env python
import os
import sys
import urllib
import urllib2
from xml.etree import ElementTree

from optparse import OptionParser

BASE_URL = "http://www.ncrna.org/frnadb/api/"
DEFAULT_NCRNA = ["tRNA","snRNA","snoRNA","rRNA","miRNA","scRNA","RNase_MRP_RNA","RNase_P_RNA","SRP_RNA","telomerase_RNA","vault_RNA","Y_RNA"]

def main(base_url, so_name, organism="Homo sapiens"):
    
    if base_url is None:
        base_url=BASE_URL
    
    if len(so_name) == 0:
        so_name = DEFAULT_NCRNA
            
    for so in so_name:
        entry_ids = {}
        sys.stderr.write("Searching for %s in %s\n" % (so,organism))
        
        for id in search(base_url,organism,[so]):
            entry_ids[id] = 1
        #for id in entry_ids.keys():
        #    print id
        total = len(entry_ids.keys())
        processed = 0
        fetched = 0
        discarded = 0
        sys.stderr.write("Found %d entries\n" % total)
        for entry_id in entry_ids.keys():
            for entry in fetch_entry(base_url, entry_id):
                s = entry['sequence'].upper().replace("U","T")
                ok = True
                for c in s:
                    if c == "A" or c == "C" or c == "G" or c == "T" or c == "N": continue
                    ok = False
                    discarded += 1
                    break
                if ok:
                    print ">%s\n%s" % ("|".join([entry['id'],entry['ontology'],entry['organism'],entry['description']]),entry['sequence'])
                    fetched += 1
                processed += 1
            if processed%100 == 0:
                sys.stderr.write("%d of %d entries processed\n" % (processed,total))
        sys.stderr.write("%d %s sequences discarded because of illegal nucleotides\n" % (discarded,so))
        sys.stderr.write("%d %s sequences fetched for %s\n" % (fetched,so,organism))
        
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
    parser.add_option("-o", "--organism", dest="organism", default=None)
    parser.add_option("-s", "--so_name", dest="so_name", action="append", default=[])
    parser.add_option("-u", "--base_url", dest="base_url", default=None)
    options, args = parser.parse_args()
    
    main(options.base_url, options.so_name, options.organism)

