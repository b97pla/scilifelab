#!/bin/bash
#SBATCH -p devel
#SBATCH -N 1
#SBATCH -t 00:00:30
#SBATCH -J pmd5deep
#SBATCH -A <project_id>

function usage {
    echo "Usage: $0 files_to_checksum"
}

# ignore directory errors redirecting to null
find | parallel md5deep {} 2> /dev/null

#usage
