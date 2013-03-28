#! /bin/bash

ROOT="/ssUppnexZone/proj"
CLEAN="0"
DRYRUN="0"
PROJECT=""
FILE="*.tar.bz2"

# Parse optional command line arguments
while getopts ":hdnf:" opt; do
  case $opt in
    d)
      CLEAN="1"
      ;;
    n)
      DRYRUN="1"
      ;;
    f)
      FILE=${OPTARG}
      ;;
    h)
      echo $"
Usage: $0 [-d -n -f FILE] PROJECT

    -d          After archiving, remove source file
    -n          Dry run, do not run any actual commands
    -f FILE     Archive the specified file, default is to find files in need of archiving
       PROJECT  Store under the specified UPPNEX project. Required
" >&2
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

PROJECT=${@:$OPTIND:1}
if [ -z ${PROJECT} ]
then
  echo "*** ERROR *** An UPPNEX project must be specified"
  exit 1
fi


# Switch to the project account, will also make sure that the iRODS module has been loaded
icd ${ROOT}/${PROJECT}
if [ $? -ne 0 ]
then
  echo "*** ERROR *** Could not icd to ${ROOT}/${PROJECT}"
  exit 1
fi

# Find files that have not been modified in the last day and process them
find ./ -maxdepth 1 -name "${FILE}" -mtime +1 |while read RUN
do
  
  RUN=`basename ${RUN}`
  DSTPATH="${ROOT}/${PROJECT}/${RUN}"
  
  # If the run already seems to have been archived, we skip it
  if [ -e ${RUN}.irods.log ]
  then
    echo `date -u +"%xT%TZ"`"   ${RUN} has already been archived. Skipping."
    continue
  fi
  
  LOGFILE=${RUN}.archive.log
  
  # Check if the file already exists on the destination, in which case we will warn and skip it
  CMD="ils ${DSTPATH}"
  ${CMD} >& /dev/null
  if [ $? -eq 0 ]
  then
    echo `date -u +"%xT%TZ"`"   ERROR: ${DSTPATH} already exists on destination, you will need to remove it with 'irm' before archiving" |tee -a ${LOGFILE}
    continue
  fi
  
  # We require a pre-computed md5-sum (computed on the original source)
  # This is just to double-check that the initial transfer was OK 
  if [ ! -e ${RUN}.md5 ]
  then
    echo `date -u +"%xT%TZ"`"   Missing pre-computed md5 sum, will NOT archive ${RUN}" |tee -a ${LOGFILE}
    continue
  fi
  
  # Checking that the pre-computed md5 matches the file
  echo `date -u +"%xT%TZ"`"   Verifying md5sum of ${RUN}" |tee -a ${LOGFILE}
  CMD="md5sum -c ${RUN}.md5"
  echo ${CMD} |tee -a ${LOGFILE}
  if [ ${DRYRUN} == "0" ]
  then
    $CMD
    if [ $? -ne 0 ]
    then
      echo `date -u +"%xT%TZ"`"   md5sum verification failed, will NOT archive ${RUN}" |tee -a ${LOGFILE}
      continue
    fi
  fi
  
  # Archiving file through iRODS
  echo `date -u +"%xT%TZ"`"   Archiving ${RUN}" |tee -a ${LOGFILE}
  CMD="iput -K -P ${RUN} ${DSTPATH}"
  echo ${CMD} |tee -a ${LOGFILE}
  if [ ${DRYRUN} == "0" ]
  then
    $CMD >& ${RUN}.irods.log
    if [ $? -ne 0 ]
    then
      echo `date -u +"%xT%TZ"`"   ERROR: Archiving of ${RUN} FAILED" |tee -a ${LOGFILE}
      mv ${RUN}.irods.log ${RUN}.failed.log
      continue
    fi
  fi
  
  # Verifying that the checksum in Swestore matches
  echo `date -u +"%xT%TZ"`"   Verifying checksum in Swestore" |tee -a ${LOGFILE}
  if [ ${DRYRUN} == "0" ]
  then
    s=`ichksum ${DSTPATH}`
    s=`echo $s |cut -f 2 -d ' '`
    u=`cut -f 1 -d ' ' ${RUN}.md5`
    if [ "$s" == "$u" ]
    then
      echo `date -u +"%xT%TZ"`"   Swestore checksum and uppmax checksum match for ${RUN}. Archiving OK" |tee -a ${LOGFILE}
    else
      echo `date -u +"%xT%TZ"`"   Swestore checksum did not match. Archiving of ${RUN} FAILED" |tee -a ${LOGFILE}
      mv ${RUN}.irods.log ${RUN}.failed.log
      continue
    fi
  else
    echo `date -u +"%xT%TZ"`"   DRYRUN: No md5 sums to verify" |tee -a ${LOGFILE}
  fi 
  
  # As an extra precaution, verify the integrity of the file using ssverify.sh
  echo `date -u +"%xT%TZ"`"   Verifying upload through ssverify.sh" |tee -a ${LOGFILE}
  CMD="ssverify.sh ${RUN} ${DSTPATH}"
  echo $CMD |tee -a ${LOGFILE}
  if [ ${DRYRUN} == "0" ]
  then
    $CMD
    if [ $? -ne 0 ]
    then
      echo `date -u +"%xT%TZ"`"   ERROR: ssverify.sh for ${RUN} FAILED" |tee -a ${LOGFILE}
      mv ${RUN}.irods.log ${RUN}.failed.log
      continue
    fi
    echo `date -u +"%xT%TZ"`"   ssverify.sh for ${RUN} PASSED. Archiving OK" |tee -a ${LOGFILE}
  fi
  
  # Remove source file if everything seems OK
  if [ ${CLEAN} == "1" ]
  then
    echo `date -u +"%xT%TZ"`"   Removing source file ${RUN}" |tee -a ${LOGFILE}
    CMD="rm -f ${RUN}"
    echo $CMD |tee -a ${LOGFILE}
    if [ ${DRYRUN} == "0" ]
    then
      $CMD
    fi
  fi
  
done
