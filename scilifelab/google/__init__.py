"""Code for writing to google docs
"""
import os


def _from_unicode(unistr, encoding='utf-8'):
    """Encode a unicode string, by default into utf-8"""
    if isinstance(unistr, unicode):
        unistr = unistr.encode(encoding)
    return unistr


def _to_unicode(str, encoding='utf-8'):
    """Decode a string into unicode"""
    if isinstance(str, basestring):
        if not isinstance(str, unicode):
            str = unicode(str, encoding)
    return str


def get_credentials(encoded_credentials_file):
    """Get the encoded credentials from the supplied file"""

    encoded_credentials = None
    if encoded_credentials_file is not None and os.path.exists(encoded_credentials_file):
        with open(encoded_credentials_file) as fh:
            encoded_credentials = fh.read().strip()
    return encoded_credentials
