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

    HLA_index = mid_HLA_index[hla_relative_index]

    Aucs = []

    for split_i in range(10):

        # load the HLA_agnostic model pvalues from corresponding split
        all_tcr_pvals = pd.read_csv(f"../HLA_agnostic_model_split/results/pvalues_split{split_i}.csv", header=0)
        all_tcr_pvals = np.array(all_tcr_pvals.pval)
        base_tcr_pool,_ = filter_p_vals(0.001, all_tcr_pvals)
        
        # load the HLA_specific model pvalues from corresponding split
        # first tell whether the current split is a valid split for the given HLA 
        indicator_i = split_mat[hla_relative_index][split_i]
        if indicator_i == 0:
            auc_i = np.nan
        else:
            test_index = test_index_mat[split_i]

            pval_dir = f"../HLA_specific_model_split/results/pvals_hla_index_{HLA_index}"
            pval_filename= os.path.join(pval_dir, f"pvalues_hla_index_{HLA_index}_{split_i}.csv")
            pvals = pd.read_csv(pval_filename, header=0)
            pvals = np.array(pvals.pval)
            important_index, _ = filter_p_vals(0.001, pvals)

            combined_TCR_index = np.concatenate([base_tcr_pool, important_index]) 
            combined_TCR_index = np.unique(combined_TCR_index) 
            
            pred_probs = make_pred(TCR, combined_TCR_index, test_index)

            auc_i = constrained_pred(HLA, HLA_index, CMV, test_index, pred_probs)

        Aucs += [auc_i]

    auc_filename = os.path.join(res_path, f"combined_split_aucs_hla_index_{HLA_index}.csv")
    
    df_auc = pd.DataFrame(list(zip(list(range(10)), Aucs)), 
                          columns=['split_id', 'AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")




