# Samples - either a space-separated list or a file name containing
# one sample name per row
samples=${samples}

# Input directory in which to look for files - make sure this points
# to the directory under which all input data can be found
indir=${indir}

# Target definition files, on interval_list format
# http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html
BAIT_INTERVALS_FILE=${baits_file}
TARGET_INTERVALS_FILE=${targets_file}

# Point to intervallist-file or define regiones for ROI. Mostly identical to target file.
TARGET_REGION=${target_region}

# Log files
${'LOGFILE={}'.format(output) if output else '# LOGFILE'}
${'ERRFILE={}'.format(error) if error else  '# ERRFILE'}
