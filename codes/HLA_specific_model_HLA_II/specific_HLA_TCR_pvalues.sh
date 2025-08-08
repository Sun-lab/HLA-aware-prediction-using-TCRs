#!/bin/bash

mkdir specific_HLA_TCR_pvalues_log

for hla_id in {0..134}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o specific_HLA_TCR_pvalues_log/specific_HLA_TCR_pvalues_${hla_id}.log \
         -p campus-new --job-name=${hla_id}_ii \
         specific_HLA_TCR_pvalues.py --hla_id ${hla_id}
done

