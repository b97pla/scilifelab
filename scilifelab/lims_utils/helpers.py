#!/usr/bin/env python

from datetime import date

def comp_dates(a, b):
    """Dates in isoformat. Is a < b?"""
    a = date(*map(int, a.split('-') ))
    b = date(*map(int, b.split('-') ))
    delta = a - b
    if delta.days < 0:
        return True
    else:
        return False

def delete_Nones(dict):
    "Deletes None type items from dict."
    new_dict = {}
    for key, val in dict.items():
        if val:
            if not val=='null':
                if not (val=='2000-10-10' or val=='3000-10-10'):
                    new_dict[key] = val
    if new_dict != {}:
        return new_dict
    else:
        return None
