#!/bin/bash

mkdir combined_HLA_log

sbatch -t 03-10:00:00 -c 1 \
        --mem=20000 --constraint=gizmok -o combined_HLA_log/combined_HLA.log \
        -p campus-new --job-name=combined combined_HLA.py


