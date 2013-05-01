#! /bin/bash

## Package a run folder to prepare it for shipping to swestore

set -e

ROOT=/srv/illumina
CLEAN="0"
SCRIPT=${HOME}/bin/convertBclToFastq.sh
DRYRUN="0"
REMOTE=""
MD5STATUS="${HOME}/swestore_transfers.log"

# Parse optional command line arguments
while getopts ":hdni:" opt; do
  case $opt in
    d)
      CLEAN="1"
      ;;
    n)
      DRYRUN="1"
      ;;
    i)
      ROOT=${OPTARG}
      ;;
    h)
      echo $"
Usage: $0 [-d -n -i ROOTDIR] RUN

    -d             After tarballing, remove original root directory
    -n             Dry run, do not run any actual commands
    -i ROOTDIR     Use the supplied path as the root directory for runs, default is ${ROOT}
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

if [ ${DRYRUN} -eq "1" ]
then
  echo "##### DRY-RUN #####"
fi

RUN=${@:$OPTIND:1}
if [ -z ${RUN} ]
then
  echo "*** ERROR *** A run name must be specified"
  exit 1
fi

REMOTEHOST=${@:$OPTIND+1:1}
REMOTEPATH=${@:$OPTIND+2:1}

if [ ! -z ${REMOTEHOST} ] && [ -z ${REMOTEPATH} ]
then
  echo "*** ERROR *** A remote path must be specified"
  exit 1
fi

STARTPATH=`pwd`
cd ${ROOT}

RUN=`basename ${RUN}`

if [ ! -e ${RUN} ]
then
  echo "No such run: ${RUN}"
  exit 1
fi

TARBALL=${RUN}.tar.bz2
HASHFILE=${TARBALL}.md5

echo `date`$'\t'"Copying convert script ${SCRIPT} to ${RUN}"
CMD="cp ${SCRIPT} ${RUN}/"
if [ ${DRYRUN} -eq "0" ]
then
  ${CMD} 
else
  echo "${CMD}"
fi

echo `date`$'\t'"Tarballing ${RUN} to ${TARBALL}"
CMD="tar --use-compress-program=pbzip2 -cf ${TARBALL} ${RUN}"
if [ ${DRYRUN} -eq "0" ]
then
  ${CMD} 
else
  echo "${CMD}"
fi

echo `date`$'\t'"Calculating md5sum of ${TARBALL} to ${HASHFILE}"
CMD="md5sum ${TARBALL}"
if [ ${DRYRUN} -eq "0" ]
then
  ${CMD} > ${HASHFILE}
else
  echo "${CMD} > ${HASHFILE}"
fi


if [ ${CLEAN} -eq "1" ]
then
  echo `date`$'\t'"Removing original run folder ${RUN}"
  CMD="rm -rf ${RUN}"
  if [ ${DRYRUN} -eq "0" ]
  then
    ${CMD} 
  else
    echo "${CMD}"
  fi
fi

if [ ! -z ${REMOTEHOST} ]
then
  REMOTE="${REMOTEHOST}:${REMOTEPATH}/"
  echo `date`$'\t'"Uploading ${TARBALL} to ${REMOTE}"
  HASHNAME=`basename ${HASHFILE}` 
  CMD1="rsync -cav ${HASHFILE} ${REMOTE}"
  CMD2="rsync -cav ${TARBALL} ${REMOTE}"
  CMD3="ssh ${REMOTEHOST} 'cd ${REMOTEPATH} && md5sum -c ${HASHNAME}'"
  CMD3=$"ssh ${REMOTEHOST} <<'ENDSSH' cd ${REMOTEPATH} && md5sum -c ${HASHNAME} ENDSSH"
  if [ ${DRYRUN} -eq "0" ]
  then
    ${CMD1}
    ${CMD2}
    ${CMD3} >> ${MD5STATUS}  
  else
    echo "${CMD1} && ${CMD2} && ${CMD3}"
  fi
fi

cd ${STARTPATH}

echo `date`$'\t'"All done!"

if [ ${DRYRUN} -eq "1" ]
then
  echo "##### DRY-RUN #####"
fi
