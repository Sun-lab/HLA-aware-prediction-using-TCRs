#!/bin/bash

mkdir HLA_agnostic_log

sbatch -t 03-10:00:00 -c 1 \
        --mem=20000 --constraint=gizmok -o HLA_agnostic_log/HLA_agnostic_base.log \
        -p campus-new --job-name=base HLA_agnostic_base.py


