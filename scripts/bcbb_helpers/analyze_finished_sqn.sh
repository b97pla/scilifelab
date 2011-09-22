#!/bin/sh
source $HOME/.bashrc
python $HOME/opt/bcbb/nextgen/scripts/analyze_finished_sqn.py $HOME/config/universe_wsgi.ini $HOME/config/post_process.yaml
