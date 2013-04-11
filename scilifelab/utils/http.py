"""Utilities for working with urls"""
import httplib
import urlparse

## From http://pythonadventures.wordpress.com/2010/10/17/check-if-url-exists/
def get_server_status_code(url):
    """
    Download just the header of a URL and
    return the server's status code.
    """
    pcs = urlparse.urlparse(url)
    try:
        conn = httplib.HTTPConnection(host=pcs.hostname, port=pcs.port)
        conn.request('HEAD', pcs.path)
        return conn.getresponse().status
    except StandardError:
        return None

def check_url(url):
    """
    Check if a URL exists without downloading the whole file.
    We only check the URL header.
    """
    good_codes = [httplib.OK, httplib.FOUND, httplib.MOVED_PERMANENTLY]
    return get_server_status_code(url) in good_codes
