#!/bin/bash

mkdir union_all_HLAs_log

sbatch -t 03-10:00:00 -c 1 \
        --mem=20000 --constraint=gizmok -o union_all_HLAs_log/union_all_HLAs.log \
        -p campus-new --job-name=union union_all_HLAs.py 


