#!/bin/bash

mkdir split_specific_model_log

for hla_relative_index in {0..20}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o split_specific_model_log/split_specific_model_hla_relative_index_${hla_relative_index}.log \
         -p campus-new --job-name=hla_${hla_relative_index} \
         split_specific_model.py --hla_relative_index ${hla_relative_index}
done

