#!/bin/bash

mkdir weighted_pred_log

for hla_id in {0..84}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o weighted_pred_log/weighted_pred_${hla_id}.log \
         -p campus-new --job-name=${hla_id}_i \
         weighted_pred.py --hla_id ${hla_id}
done

