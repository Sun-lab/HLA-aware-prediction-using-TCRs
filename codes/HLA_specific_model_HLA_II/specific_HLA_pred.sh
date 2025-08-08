#!/bin/bash

mkdir specific_HLA_pred_log

for hla_id in {0..134}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o specific_HLA_pred_log/specific_HLA_pred_${hla_id}.log \
         -p campus-new --job-name=${hla_id}_i \
         specific_HLA_pred.py --hla_id ${hla_id}
done

