#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import multiprocessing
import numpy as np
from scipy import io
import os
import pandas as pd
from scipy.stats import fisher_exact
from sklearn import metrics
from datetime import datetime
import argparse

def read_in_data(data_path):
    exception_index = [60, 168,   7,  88, 131]
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    split_mat = np.array(pd.read_csv(os.path.join(data_path, "split_mat.csv"), header=None, dtype=int))
    mid_HLA_index = np.array(pd.read_csv(os.path.join(data_path, "mid_HLA_index.csv"), header=None, dtype=int)[0])
    train_index_mat = np.array(pd.read_csv(os.path.join(data_path, "split_train_index.csv"), header=None, dtype=int))
    test_index_mat = np.array(pd.read_csv(os.path.join(data_path, "split_test_index.csv"), header=None, dtype=int))

    new_HLA = np.vstack((HLA, HLA[exception_index]))

    return  TCR, new_HLA, CMV, split_mat, mid_HLA_index, train_index_mat, test_index_mat


def get_fisher_data(training_index, TCR):
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


def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    return res


def filter_p_vals(cutoff, p_vals):
    p_val_res = np.zeros(len(p_vals))
    index_res = np.zeros(len(p_vals), int)
    res_index = 0
    for i in range(len(p_vals)):
        if p_vals[i] <= cutoff:
            p_val_res[res_index] = p_vals[i]
            index_res[res_index] = int(i)
            res_index += 1
    
    return index_res[0:res_index],  p_val_res[0: res_index], res_index



def make_pred(TCR, important_index, test_individuals):
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = [ testing_TCR[i] for i in important_index ]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return prediction_probs



def constrained_pred(HLA, target_HLA, CMV, test_indexes, pred_probs, split_info):
    AUC = []
    for posi_index, HLA_i in enumerate(target_HLA):
        if split_info[posi_index] == 0:
            AUC.append(np.nan)
            continue
        position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
        HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
        position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
        if len(HLA_i_test_index)<=1:
            Auc_i = np.nan
        elif len(set(CMV[HLA_i_test_index]))==1:
            Auc_i = np.nan
        else:
            Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
        AUC.append(Auc_i)
    return AUC


def work_wrapper(TCR, HLA, CMV, train_index_mat, test_index_mat, type1_error, split_index, split_mat, mid_HLA_index):

    cardi = len(TCR)
    train_index = train_index_mat[split_index]
    test_index = test_index_mat[split_index]

    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_index, TCR)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
    
    significant_TCR_index, signi_TCR_pvals, sign_count = filter_p_vals(type1_error, p_vals)

    pred_probs = make_pred(TCR, significant_TCR_index, test_index )

    Auc = constrained_pred(HLA, mid_HLA_index, CMV, test_index, pred_probs, split_mat[:, split_index])

    return Auc, p_vals



# parellize the 10 different splits, this is because for each split, we only need to train the model once
# within each split's subroutine, we loop over 21 HLAs in constrained predcition

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run HLA-agnostic model on all training samples, under each given one of 10 splits')
    parser.add_argument('--split_id', type=int, help='the index out of the ten splits')
    input_args = parser.parse_args()
    split_id = input_args.split_id

    print(split_id)

    data_path = "../intermediate_files"
    res_path =  "./results"
    os.makedirs(res_path, exist_ok=True)

    auc_file_name = 'aucs_split'+str(split_id)+".csv"
    pvalue_file_name = "pvalues_split"+str(split_id)+".csv"

    type1_error = 0.001
    TCR, HLA, CMV, split_mat, mid_HLA_index, train_index_mat, test_index_mat = read_in_data(data_path)

    AUC_values, np_pvals = work_wrapper(TCR, HLA, CMV, train_index_mat, test_index_mat, type1_error, split_id, split_mat, mid_HLA_index)

    df_pvalues = pd.DataFrame(np_pvals.tolist(), 
                              columns=['pval'])
    df_pvalues.to_csv(os.path.join(res_path, pvalue_file_name), 
                      index=False, na_rep="NA")

    df_auc = pd.DataFrame(list(zip(mid_HLA_index.tolist(), 
                                   AUC_values)), 
                          columns=['HLA_index', 'AUC'])
    df_auc.to_csv(os.path.join(res_path, auc_file_name), 
                  index=False, na_rep="NA")



