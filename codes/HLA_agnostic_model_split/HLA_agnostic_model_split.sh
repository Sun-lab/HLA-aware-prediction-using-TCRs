#!/bin/bash

mkdir HLA_agnostic_model_split_log

for split_id in {0..9}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o HLA_agnostic_model_split_log/HLA_agnostic_model_split_${split_id}.log \
         -p campus-new --job-name=${split_id}_split \
         HLA_agnostic_model_split.py --split_id ${split_id}
done

