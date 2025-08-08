#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import io
import os
import pandas as pd
from sklearn import metrics




def make_pred(TCR, important_index, test_individuals):
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = testing_TCR[important_index]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return prediction_probs


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
    HLA_I_index = np.array(pd.read_csv(os.path.join(data_path, "HLA_I_index.csv"),header=None,dtype=int )[0])
    HLA_II_index = np.array(pd.read_csv(os.path.join(data_path, "HLA_II_index.csv"),header=None,dtype=int )[0])    
    return TCR, new_HLA, CMV, train_index, test_index, HLA_I_index, HLA_II_index


def constrained_pred(HLA, HLA_i, CMV, test_indexes, pred_probs):
    test_indexes = np.array(test_indexes)
    position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
    HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
    position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
    if len(HLA_i_test_index)<=1:
        Auc_i = np.nan
    elif len(set(CMV[HLA_i_test_index]))==1:
        Auc_i = np.nan
    else:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
    return Auc_i, CMV[ HLA_i_test_index ], pred_probs[position_index_pred]



if __name__ == "__main__":

    data_path = "../intermediate_files/"
    res_path ="./results" 
    roc_path ="./results/roc_related"
    os.makedirs(res_path, exist_ok=True)
    os.makedirs(roc_path, exist_ok=True)
    

    TCR, HLA, CMV, train_index, test_index, HLA_I_index, HLA_II_index = read_in_data("../intermediate_files")
    all_tcr_pvals = pd.read_csv("../HLA_agnostic_model/results/pvalues.csv", header=0)
    all_tcr_pvals = np.array(all_tcr_pvals.pval)
    base_tcr_pool,_ = filter_p_vals(0.001, all_tcr_pvals)
    all_hlas = np.arange(220)
    Aucs = []
    # this for loop reads in all the HLA pvals and create the HLA index <-> significant TCR index map
    for i in all_hlas:

        if i in HLA_I_index:
            pval_dir = "../HLA_specific_model_HLA_I/results/pvals"
        else:
            pval_dir = "../HLA_specific_model_HLA_II/results/pvals"   

        pval_filename= os.path.join(pval_dir, f"pvalues_hla_index_{i}.csv")
        pvals = pd.read_csv(pval_filename, header=0)
        pvals = np.array(pvals.pval)

        important_index, _ = filter_p_vals(0.001, pvals)

        combined_TCR_index = np.concatenate([base_tcr_pool, important_index]) 
        combined_TCR_index = np.unique(combined_TCR_index) 

        pred_probs = make_pred(TCR, combined_TCR_index, test_index)

        auc_i, true_labels, pred_numeric = constrained_pred(HLA, i, CMV, test_index, pred_probs)
        Aucs.append(auc_i)

        if auc_i==auc_i:
            roc_filename = os.path.join(roc_path, f"roc_related_{i}.csv")
            df_roc = pd.DataFrame(list(zip(true_labels.tolist(), pred_numeric.tolist())), 
                                columns=["label", "prediction"])
            df_roc.to_csv(roc_filename, index=False, na_rep="NA")

    auc_filename = os.path.join(res_path, "combined_aucs.csv")
    
    df_auc = pd.DataFrame(list(zip(list(range(HLA.shape[0])), Aucs)), 
                          columns=['HLA_index', 'AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")
    

