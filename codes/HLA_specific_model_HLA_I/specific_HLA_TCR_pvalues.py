#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
for each HLA, output the association between the given HLA and all TCRs among the training individuals
need to handle the case that some individuals 
do not have HLA typing results for some HLA locations
"""


import multiprocessing
import numpy as np
from scipy import io
import pandas as pd
import os
from scipy.stats import fisher_exact
from sklearn import metrics
from datetime import datetime
import argparse


def get_fisher_data(training_index, TCR, HLA_row):
    cardi = len(TCR)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    res_index = 0
    
    for i in range(cardi):
        for j in training_index:
            if TCR[i,j] == 1 and HLA_row[j] == 1 :
                oneone_vec[res_index] += 1
                continue
            if TCR[i,j] == 1 and HLA_row[j] == 0 :
                onezero_vec[res_index] += 1
                continue
            if TCR[i,j] == 0 and HLA_row[j] == 0 :
                zerozero_vec[res_index] += 1
                continue
            if TCR[i,j] == 0 and HLA_row[j] == 1 :
                zeroone_vec[res_index] += 1
                continue
        res_index += 1
    return oneone_vec, onezero_vec, zeroone_vec, zerozero_vec


def read_in_data(data_path, intm_path):
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_matrix = pd.read_csv( os.path.join(intm_path, "complete_HLA_NAs_kept.csv") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_matrix.to_numpy()
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    AA_index = np.array(pd.read_csv(os.path.join(data_path, "AA_HLA_index.csv"),header=None,dtype=int )[0])

    
    return TCR, HLA, train_index, test_index, AA_index



def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    return res
    

def work_wrapper(TCR, HLA, AA_index, train_index, test_indexes, pid):
    cardi = len(TCR)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    target_HLA = AA_index[pid]
    HLA_row = HLA[target_HLA,]
    HLA_known_individuals = np.where(~np.isnan(HLA_row))[0]
    train_indexes = np.intersect1d(train_index, HLA_known_individuals)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR, HLA_row)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)

    return p_vals

# no clustering, for each HLA, we use the subset of individuals from the training data 
# who have this HLA to compute the association pvalue between each TCR and CMV status

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='association pvalues between TCR and CMV under each HLA')
    parser.add_argument('--hla_id', type=int, help='the index of HLA (among HLA-Is) to get TCR and CMV association p-values for')
    input_args = parser.parse_args()
    hla_id = input_args.hla_id

    print(hla_id)

    data_path ="./data_files"
    intm_path ="../intermediate_files"
    res_path ="./results/TCR_HLA_pvalues"
    os.makedirs(res_path, exist_ok=True)

    res_file_name = str(hla_id)+"_TCR_HLA_pvalues.csv"

    #data_path = "/home/hyo/mystuff/hutch_project/data"
    TCR,HLA,train_index,test_index,AA_index = read_in_data(data_path, intm_path)
    
    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%H:%M:%S")
    print(current_time)
    print("start computing p-values")

    np_pvalues = work_wrapper(TCR, HLA, AA_index, train_index, test_index, hla_id)

    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%H:%M:%S")
    print(current_time)
    print("finished computing p-values")

    df_pvalues = pd.DataFrame(np_pvalues.tolist(), 
                              columns=['pval'])
    
    df_pvalues.to_csv(os.path.join(res_path, res_file_name), 
                            index=False)

