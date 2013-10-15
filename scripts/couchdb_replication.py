#!/usr/bin/env python

import couchdb
import json
import argparse
import logbook
import sys

from couchdb import PreconditionFailed

#Set up logging
l = logbook.Logger(level=logbook.INFO)
h = logbook.StreamHandler(sys.stdout, level=logbook.INFO)


def _get_databases_info(source, destination):
    """Returns a tuple containing a python representation of source and destination
    couchDB instances. It also returns a list of the databases in both instances.
    """
    s_couch = couchdb.Server(source)
    d_couch = couchdb.Server(destination)
    _, _, s_dbs = s_couch.resource.get_json('_all_dbs')
    _, _, d_dbs = d_couch.resource.get_json('_all_dbs')
    return s_couch, d_couch, s_dbs, d_dbs



def _setup_continuous(source, destination):
    """Set up a continuous replication of all databases in source to destination.
    """
    s_couch, d_couch, s_dbs, d_dbs = _get_databases_info(source, destination)

    l.info("Databases in the source CouchDB instance: {}".format(', '.join(s_dbs)))
    l.info("Databases in the destination CouchDB instance: {}".format(', '.join(d_dbs)))

    #We don't want to replicate the replicator DB
    try:
        s_dbs.remove('_replicator')
        d_dbs.remove('_replicator')
    except ValueError:
        pass

    #For each DB in the source CouchDB instance, create a replication document
    #and get its _security object to put it in the destination database
    for db in s_dbs:
        _, _, security = s_couch[db].resource.get_json('_security')
        doc = {
                'name': '{}_rep'.format(db),
                'source': '{}/{}/'.format(source, db),
                'target': '{}/{}/'.format(destination, db),
                'continuous': True
        }
        s_rep = s_couch['_replicator']

        #Create the DB in the destination if not present
        try:
            d_couch.create(db)
            l.info("Created {} database in destination".format(db))
        except PreconditionFailed:
            l.info("Database {} already existing in the destination, not creating it".format(db))

        #Put the replicator document in source and set security object in destination
        l.info("Putting replicator document in _replicator database of source")
        s_rep.create(doc)
        l.info("Copying security object to {} database in destination".format(db))
        d_couch[db].resource.put('_security', security)

        l.info("DONE!")



if __name__ == "__main__":

    DESCRIPTION = """Set up complete one-way replication for CouchDB.

    Use this script if you want to configure a stage database that will have the
    exact same content of your production database.

    To do so, the script creates a replication document for each database in the
    source CouchDB instance that replicates such database (in continuous mode)
    to the destination database.

    Security object (permissions per database), are put to the destination databases.
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('action', type=str, help = "Action to perform, either \
            configure continuous replication (continuous) or punctual clone (clone)")
    parser.add_argument('source', type=str, help = "Source CouchDB instance, \
            with the credentials included in the URL. I.E: http://admin:passw@source_db:5984")
    parser.add_argument('destination', type=str, help = "Destination CouchDB instance, \
            with the credentials included in the URL. I.E: http://admin:passw@destination_db:5984")

    args = parser.parse_args()
    source = args.source
    destination = args.destination
    action = args.action

    with h.applicationbound():
        actions = ['continuous', 'clone']
        if action not in actions:
            raise ValueError("Action not recognised, please choose between %s" % \
                    ', '.join(actions))
        if action == "continuous":
            _setup_continuous(source, destination)
        else:
            _clone(source, destination)

