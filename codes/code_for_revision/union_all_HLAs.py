#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import io
import os
import pandas as pd
from sklearn import metrics
from collections import defaultdict
from collections import Counter



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


def predict_test_individual(test_index, sign_dict, HLA, TCR):
    pred_probs = []
    
    test_index_list = test_index.tolist()
    for ind_index in test_index_list:
        existing_hlas = np.where(HLA[:,ind_index]==1)[0]
        ind_sign_union = []
        for i in existing_hlas:
            ind_sign_union += sign_dict[i]
        ind_sign_union = list(set(ind_sign_union))
        ind_sign_union.sort()
        ind_sign_union = np.array(ind_sign_union)
        
        if len(ind_sign_union) == 0:
            pred_probs.append(np.float64(0))
        else:
            prop = make_pred(TCR, ind_sign_union, ind_index)
            pred_probs.append(prop)
        
    return np.array(pred_probs)


if __name__ == "__main__":

    data_path = "../intermediate_files/"
    res_path ="./results" 
    os.makedirs(res_path, exist_ok=True)

    TCR, HLA, CMV, train_index, test_index, HLA_I_index, HLA_II_index = read_in_data("../intermediate_files")

    # build a dictionary with HLA index as key and the list of corresponding significant TCRs as value
    sign_dict = defaultdict(list)

    for i in range(220):
        if i in HLA_I_index:
            pval_dir = "../HLA_specific_model_HLA_I/results/pvals"
        else:
            pval_dir = "../HLA_specific_model_HLA_II/results/pvals"   

        pval_filename= os.path.join(pval_dir, f"pvalues_hla_index_{i}.csv")
        pvals = pd.read_csv(pval_filename, header=0)
        pvals = np.array(pvals.pval)

        important_index, _ = filter_p_vals(0.001, pvals)
        sign_dict[i] = important_index.tolist()

    print(len(sign_dict))
    print(Counter([len(sign_dict[x])>0 for x in range(220)]))

    pred_probs = predict_test_individual(test_index, sign_dict, HLA, TCR)

    Auc = metrics.roc_auc_score( CMV[test_index ], pred_probs )

    auc_filename = os.path.join(res_path, "union_all_HLAs_auc.csv")
    
    df_auc = pd.DataFrame([Auc], 
                          columns=['AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")
    

