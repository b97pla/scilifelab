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


def get_credentials(config):
    """Get the encoded credentials specified in the post process configuration file"""

    encoded_credentials = None
    gdocs = config.get("gdocs_upload", {})
    encoded_credentials_file = gdocs.get("gdocs_credentials", None)
    if encoded_credentials_file is not None and os.path.exists(encoded_credentials_file):
        with open(encoded_credentials_file) as fh:
            encoded_credentials = fh.read().strip()
    return encoded_credentials
