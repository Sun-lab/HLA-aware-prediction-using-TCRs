#!/bin/bash

mkdir split_knn_log

for hla_id in {0..20}; do
  for split_id in {0..9}; do
  sbatch -t 03-10:00:00 -c 1 \
         --mem=20000 --constraint=gizmok -o split_knn_log/split_knn_${hla_id}_${split_id}.log \
         -p campus-new --job-name=${hla_id}_${split_id} \
         split_knn.py --hla_id ${hla_id} --split_id ${split_id}
  done
done

