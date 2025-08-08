#!/bin/bash

mkdir combined_HLA_split_log

for hla_relative_index in {0..20}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o combined_HLA_split_log/combined_HLA_split_hla_relative_index_${hla_relative_index}.log \
         -p campus-new --job-name=hla_${hla_relative_index} \
         combined_HLA_split.py --hla_relative_index ${hla_relative_index}
done

