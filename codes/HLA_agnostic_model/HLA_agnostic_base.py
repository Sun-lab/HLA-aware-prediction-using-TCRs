#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from scipy import io
import os
import pandas as pd
from scipy.stats import fisher_exact
from sklearn import metrics


def read_in_data(data_path):
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    return TCR, CMV, train_index, test_index


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
    
    return index_res[0:res_index],  p_val_res[0: res_index]



def make_pred(TCR, important_index, test_individuals):
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = [ testing_TCR[i] for i in important_index ]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return prediction_probs


def work_wrapper(TCR, CMV, train_index, test_index, cutoffs):
    cardi = len(TCR)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_index, TCR, CMV)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
    Aucs = []
    for cutoff in cutoffs:
        significant_TCR_index, signi_TCR_pvals = filter_p_vals(cutoff, p_vals)
        pred_probs = make_pred(TCR, significant_TCR_index, test_index )
        Auc = metrics.roc_auc_score( CMV[test_index ], pred_probs )
        Aucs += [Auc]
    return Aucs, p_vals



if __name__ == "__main__":
    
    data_path = "../intermediate_files/"
    res_path ="./results" 
    os.makedirs(res_path, exist_ok=True)

    cutoffs = [0.01, 0.001, 0.0001, 0.00001]
    TCR, CMV, train_index, test_index = read_in_data(data_path)
    Aucs, p_vals = work_wrapper(TCR, CMV, train_index, test_index, cutoffs)

    pval_filename = os.path.join(res_path, "pvalues.csv")

    df_pval = pd.DataFrame(p_vals.tolist(), 
                           columns=['pval'])
    df_pval.to_csv(pval_filename, index=False, na_rep="NA")


    auc_filename = os.path.join(res_path, "base_aucs.csv")

    df_auc = pd.DataFrame(list(zip(cutoffs, Aucs)), 
                          columns=['pval_cutoff', 'AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")
    
