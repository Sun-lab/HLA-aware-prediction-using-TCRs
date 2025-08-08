#!/bin/bash

mkdir specific_HLA_pred_permutation_log

for hla_id in 23; do
for chunk_id in {0..99}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o specific_HLA_pred_permutation_log/specific_HLA_pred_permutation_${hla_id}_${chunk_id}.log \
         -p campus-new --job-name=${hla_id}_${chunk_id} \
         specific_HLA_pred_permutation.py --hla_id ${hla_id} --chunk_id ${chunk_id}
done
done
