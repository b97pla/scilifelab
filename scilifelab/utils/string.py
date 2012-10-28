"""String utils module"""
import re


def replace_ascii(str):
    # Substitute swedish characters for sensible counterparts
    str = str.replace(u'\xc5','A')
    str = str.replace(u'\xe5','a')
    str = str.replace(u'\xc4','A')
    str = str.replace(u'\xe4','a')
    str = str.replace(u'\xd6','O')
    str = str.replace(u'\xf6','o') 
    return str.encode('ascii','replace')

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings of equal lengths.
    Raise ValueError if strings are of unequal length.
    """
    if len(s1) != len(s2): raise ValueError('strings of unequal length')
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))    

def strip_extensions(fn, ext=[]):
    """Strip extensions from a filename.

    :param fn: filename
    :param ext: list of extensions to strip

    :returns: stripped version of fn and extension
    """
    pattern = "|".join(".*({})$".format(x) for x in ext)
    if not pattern:
        return (fn, None)
    m = re.search(pattern, fn)
    if not m:
        return (fn, None)
    return (fn.replace(m.group(1), ""), m.group(1))
