#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
For individuals with both HLA-A*02:01 and HLA-B*07:02
'''

import multiprocessing
import numpy as np
from scipy import io
import pandas as pd
import os
from scipy.stats import fisher_exact
from sklearn import metrics
import argparse


def filter_p_vals(cutoff, p_vals):
    p_val_res = np.zeros(len(p_vals))
    index_res = np.zeros(len(p_vals), int)
    res_index = 0
    for i in range(len(p_vals)):
        if p_vals[i] <= cutoff:
            p_val_res[res_index] = p_vals[i]
            index_res[res_index] = int(i)
            res_index += 1
    
    return index_res[0:res_index],  p_val_res[0: res_index]


def make_pred(TCR, important_index, test_individuals):
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = [ testing_TCR[i] for i in important_index ]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return prediction_probs



def get_fisher_data(training_index, TCR, CMV):
    cardi = len(TCR)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    res_index = 0
    

    for i in range(cardi):
        for j in training_index:
            if TCR[i,j] == 1 and CMV[j] == 1 :
                oneone_vec[res_index] += 1
                continue
            if TCR[i,j] == 1 and CMV[j] == 0 :
                onezero_vec[res_index] += 1
                continue
            if TCR[i,j] == 0 and CMV[j] == 0 :
                zerozero_vec[res_index] += 1
                continue
            if TCR[i,j] == 0 and CMV[j] == 1 :
                zeroone_vec[res_index] += 1
                continue
        res_index += 1
    return oneone_vec, onezero_vec, zeroone_vec, zerozero_vec


def read_in_data(data_path):
    exception_index = [60, 168,   7,  88, 131]

    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    new_HLA = np.vstack((HLA, HLA[exception_index]))

    CMV = df['x']
    CMV = CMV.to_numpy()
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    HLA_names = np.array(pd.read_csv(os.path.join(data_path, "complete_HLA_rownames.csv"),header=0,dtype=str )["HLA_name"])

    
    return TCR, new_HLA, CMV, train_index, test_index, HLA_names



def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    return res


def constrained_pred(HLA, double_index, CMV, test_indexes, pred_probs):
    position_index = np.arange(len(test_indexes))
    HLA_both_status = np.sum(HLA[double_index,], axis=0)
    HLA_both_test_index = test_indexes[ HLA_both_status[test_indexes] == 2 ]
    position_index_pred = position_index[ HLA_both_status[test_indexes] == 2 ]
    if len(HLA_both_test_index)<=1:
        Auc_i = np.nan
    elif len(set(CMV[HLA_both_test_index]))==1:
        Auc_i = np.nan
    else:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_both_test_index ], pred_probs[position_index_pred] )
    return Auc_i
    

def work_wrapper(TCR, HLA, CMV, double_index, train_index, test_indexes, cutoff):
    cardi = len(TCR)
    HLA_sharing_individuals = np.where(np.sum(HLA[double_index,], axis=0)==2)[0]
    train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR, CMV)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)

    significant_TCR_index, signi_TCR_pvals = filter_p_vals(cutoff, p_vals)
    pred_probs = make_pred(TCR, significant_TCR_index, test_indexes )
    Auc_i = constrained_pred(HLA, double_index, CMV, test_indexes, pred_probs)

    return Auc_i, p_vals

# for each HLA, we build the model using all the individuals that carry this HLA.
# Recommanding parallize the HLA.
# The prediction is constrained to the test individuals that have this particular HLA.

if __name__ == "__main__":

    data_path ="../intermediate_files"
    res_path ="./results" 
    os.makedirs(res_path, exist_ok=True)

    cutoff = 0.001

    TCR,HLA,CMV,train_index,test_index,HLA_names = read_in_data(data_path)
    
    hla_index_1 = np.where(HLA_names=="HLA-A*02:01")[0].item()
    hla_index_2 = np.where(HLA_names=="HLA-B*07:02")[0].item()
    
    double_index = [hla_index_1, hla_index_2]

    Auc, p_vals = work_wrapper(TCR, HLA, CMV, double_index,
                                train_index, test_index, cutoff)


    pval_filename = os.path.join(res_path, "pvalues_double_hlas.csv")

    df_pval = pd.DataFrame(p_vals.tolist(), 
                           columns=['pval'])
    df_pval.to_csv(pval_filename, index=False, na_rep="NA")


    auc_filename = os.path.join(res_path, "auc_double_hlas.csv")

    df_auc = pd.DataFrame([Auc], 
                          columns=['AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")
    