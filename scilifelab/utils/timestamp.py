"""
Utilities for comparing time stamps
"""
import os
from datetime import datetime
import time

def utc_time():
    """Make an utc_time with appended 'Z'"""
    return str(datetime.utcnow()) + 'Z'
    
def modified_within_days(item, n_days):
    """Check if file/directory has been modified within n days.

    :param item: file/directory
    :param n_days: number of days

    :returns: boolean
    """
    if time.time()  - os.path.getmtime(item) < n_days * 3600 * 24:
        return True
    return False

