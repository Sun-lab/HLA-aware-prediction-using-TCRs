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



def constrained_pred(HLA, HLA_i, CMV, test_indexes, pred_probs):
    position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
    HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
    position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
    if len(HLA_i_test_index)<=1:
        Auc_i = np.nan
    elif len(set(CMV[HLA_i_test_index]))==1:
        Auc_i = np.nan
    else:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
    return Auc_i



def work_wrapper(TCR, HLA, CMV, train_index_mat, test_index_mat, type1_error, hla_relative_index, split_mat, mid_HLA_index, res_path):
    temp_res = []
    cardi = len(TCR)
    # HLA_index is the index of the allels in the HLA matrix
    HLA_index = mid_HLA_index[hla_relative_index]
    # prepare path to store pvalue files
    pval_path = os.path.join(res_path, f"pvals_hla_index_{HLA_index}")
    os.makedirs(pval_path, exist_ok=True)
    # prepare path to store AUC files
    auc_path = os.path.join(res_path, f"res_aucs_hla_index_{HLA_index}")
    os.makedirs(auc_path, exist_ok=True)

    for i in range(10):
        indicator_i = split_mat[hla_relative_index][i]
        train_index = train_index_mat[i]
        test_index = test_index_mat[i]
        
        if indicator_i == 0:
            temp_res.append(np.nan)
            continue
    
        HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
        train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
        oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR)
        p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
        # write out p_val file
        pval_filename = f"pvalues_hla_index_{HLA_index}_{i}.csv"
        df_pvalues = pd.DataFrame(p_vals.tolist(), 
                                  columns=['pval'])
        df_pvalues.to_csv(os.path.join(pval_path, pval_filename), 
                          index=False, na_rep="NA")

        significant_TCR_index, signi_TCR_pvals, sign_count = filter_p_vals(type1_error, p_vals)

        pred_probs = make_pred(TCR, significant_TCR_index, test_index )
        Auc_i = constrained_pred(HLA, HLA_index, CMV, test_index, pred_probs)
        temp_res.append(Auc_i)
        
    # write out the list of AUCs
    df_aucs = pd.DataFrame(temp_res, 
                           columns=["AUC"])
    auc_filename = f"aucs_hla_index_{HLA_index}.csv"
    df_aucs.to_csv(os.path.join(auc_path, auc_filename), 
                   index=False, na_rep="NA") 

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run HLA-specific model on training samples under 10 splits, under given HLA relative index')
    parser.add_argument('--hla_relative_index', type=int, help='the relative index out of the list of HLAs with frequency in the middle')
    input_args = parser.parse_args()
    hla_relative_index = input_args.hla_relative_index

    print(hla_relative_index)

    data_path = "../intermediate_files"
    res_path =  "./results"
    os.makedirs(res_path, exist_ok=True)

    type1_error = 0.001
    TCR, HLA, CMV, split_mat, mid_HLA_index, train_index_mat, test_index_mat = read_in_data(data_path)

    work_wrapper(TCR, HLA, CMV, train_index_mat,
                    test_index_mat, type1_error, hla_relative_index, split_mat, mid_HLA_index, res_path)






