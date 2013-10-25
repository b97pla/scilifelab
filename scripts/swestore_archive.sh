#! /bin/bash -l

# Function for echoing to stderr
echoerr() { echo "$@" 1>&2; }

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
  echoerr "*** ERROR *** An UPPNEX project must be specified"
  exit 1
fi


if [ ${DRYRUN} == "0" ]
then
  N=""
else
  N="-n "
fi

if [ ${CLEAN} == "0" ]
then
  D=""
else
  D="-d "
fi

DIR=`dirname ${FILE}`
FNAME=`basename ${FILE}`
# Find files and process them
find ${DIR} -maxdepth 1 -name "${FNAME}" |while read RUN
do
  CMD="bash swestore_archive_run.sh ${D}${N}${RUN} ${PROJECT}"
  echo $CMD
  $CMD
done
  