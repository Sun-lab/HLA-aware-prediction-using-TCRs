#!/usr/bin/env python3
# -*- coding: utf-8 -*-




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


def read_in_data(data_path):
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    AA_index = np.array(pd.read_csv(os.path.join(data_path, "HLA_I_index.csv"),header=None,dtype=int )[0])

    
    return TCR, HLA, CMV, train_index, test_index, AA_index



def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    return res


def constrained_pred(HLA_i, CMV, test_indexes, pred_probs):
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
    

def work_wrapper(TCR, HLA, CMV, AA_index, train_index, test_indexes, cutoffs, hla_id):
    cardi = len(TCR)
    target_HLA = AA_index[hla_id]
    HLA_sharing_individuals = np.where(HLA[target_HLA,]  > 0)[0]
    train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
    Aucs = []
    for cutoff in cutoffs:
        significant_TCR_index, signi_TCR_pvals = filter_p_vals(cutoff, p_vals)
        pred_probs = make_pred(TCR, significant_TCR_index, test_indexes )
        Auc_i = constrained_pred(target_HLA, CMV, test_indexes, pred_probs)
        Aucs += [Auc_i]
    return Aucs, p_vals

# for each HLA, we build the model using all the individuals that carry this HLA.
# Recommanding parallize the HLA.
# The prediction is constrained to the test individuals that have this particular HLA.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='HLA specific prediction')
    parser.add_argument('--hla_id', type=int, help='the relative index of HLA (among HLA-Is) to get prediction')
    input_args = parser.parse_args()
    hla_id = input_args.hla_id

    print(hla_id)

    data_path ="../intermediate_files"
    pvalue_path ="./results/pvals"
    auc_path ="./results/aucs"   
    os.makedirs(pvalue_path, exist_ok=True)
    os.makedirs(auc_path, exist_ok=True)

    cutoffs = [0.01, 0.001, 0.0001, 0.00001]

    TCR,HLA,CMV,train_index,test_index,AA_index = read_in_data(data_path)
    

    Aucs, p_vals = work_wrapper(TCR, HLA, CMV, AA_index,
                                train_index, test_index, cutoffs, hla_id)

    target_HLA = AA_index[hla_id]


    pval_filename = os.path.join(pvalue_path, f"pvalues_hla_index_{target_HLA}.csv")

    df_pval = pd.DataFrame(p_vals.tolist(), 
                           columns=['pval'])
    df_pval.to_csv(pval_filename, index=False, na_rep="NA")


    auc_filename = os.path.join(auc_path, f"aucs_{target_HLA}.csv")

    df_auc = pd.DataFrame(list(zip(cutoffs, Aucs)), 
                          columns=['pval_cutoff', 'AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")
    