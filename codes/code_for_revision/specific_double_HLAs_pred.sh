#!/bin/bash

mkdir specific_double_HLAs_pred_log


sbatch -t 03-10:00:00 -c 1 \
        --mem=20000 --constraint=gizmok -o specific_double_HLAs_pred_log/specific_double_HLAs_pred.log \
        -p campus-new --job-name=double specific_double_HLAs_pred.py 


