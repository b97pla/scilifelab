#!/bin/sh
# Re-runs the samplesheet to yaml converter and copies it on the remote analysis machine

samplesheet=$1
BCBB_PATH="$HOME/opt/bcbb/nextgen/scripts"
STORE_HOST=`grep store_host: ~/config/post_process.yaml | cut -f2 -d" "`
STORE_PATH=`grep store_dir: ~/config/post_process.yaml | cut -f2 -d" "`

python $BCBB_PATH/utils/convert_samplesheet_config.py $samplesheet

# drop file extension
filename=${samplesheet%.*}

echo "Copying $filename.yaml to $STORE_HOST..."

scp $filename.yaml "$STORE_HOST:$STORE_PATH/run_info.yaml"
rm $filename.yaml

echo "The resulting run_info.yaml file has been moved to $STORE_PATH/run_info.yaml"
