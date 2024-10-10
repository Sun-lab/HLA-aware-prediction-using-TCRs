#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from scipy import io
import os
import pandas as pd
from scipy.stats import fisher_exact
from sklearn import metrics

def read_in_data(data_path):
    exception_index = [60, 168,   7,  88, 131]
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    new_HLA = np.vstack((HLA, HLA[exception_index]))
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    return TCR, new_HLA, CMV, train_index, test_index



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
    
    return index_res[0:res_index],  p_val_res[0: res_index]



def make_pred(TCR, important_index, test_individuals):
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = [ testing_TCR[i] for i in important_index ]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return prediction_probs


def constrained_pred(target_HLA, CMV, test_indexes, pred_probs):
    test_indexes = np.array(test_indexes)
    AUC = []
    for HLA_i in target_HLA:
        position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
        HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
        position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
        try:
            Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
        except ValueError:
            Auc_i = -1
        AUC.append(Auc_i)
    return np.array(AUC)



def work_wrapper(TCR, CMV, train_index, test_index, type1_error):
    cardi = len(TCR)
    HLAs = np.arrange(220)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_index, TCR)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
    significant_TCR_index, signi_TCR_pvals = filter_p_vals(type1_error, p_vals)
    pred_probs = make_pred(TCR, significant_TCR_index, test_index )
    Auc = constrained_pred(HLAs, CMV, test_index, pred_probs)
    return Auc



if __name__ == "__main__":
    data_path = "../intermediate_files/"
    type1_error = 0.001
    TCR, CMV, train_index, test_index = read_in_data(data_path)
    res = work_wrapper(TCR, CMV, train_index, test_index, type1_error)
    np.savetxt("HLA_ignorant_specifc_res.csv", res)

