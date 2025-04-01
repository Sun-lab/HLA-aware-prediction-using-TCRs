#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 22:25:04 2024

@author: hyo
"""
import multiprocessing
import numpy as np
from scipy import io
import os
import pandas as pd
from scipy.stats import fisher_exact
from sklearn import metrics


def read_in_data(data_path):
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

    return  TCR, HLA, CMV, split_mat, mid_HLA_index, train_index_mat, test_index_mat


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
    #HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
    #training_TCR = TCR[:, HLA_sharing_individuals]
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
    try:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
    except ValueError:
        Auc_i = -1
            
    return Auc_i



def work_wrapper(TCR, HLA, CMV, train_index, test_index, type1_error, HLA_index, split_mat, res_dict):
    temp_res = np.array([])
    cardi = len(TCR)

    for indicator_i in split_mat[HLA_index]:
        if indicator_i == 0:
            np.append(temp_res, -1)
            continue
    
        oneone_vec = np.zeros(cardi)
        onezero_vec = np.zeros(cardi)
        zeroone_vec = np.zeros(cardi)
        zerozero_vec = np.zeros(cardi)
        HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
        train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
        oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR)
        p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
        significant_TCR_index, signi_TCR_pvals = filter_p_vals(type1_error, p_vals)
        #res_name1 = os.path.join(res_path, "{}TCR_pvals.txt".format(target_HLA))
        #np.savetxt(res_name1, p_vals)
        pred_probs = make_pred(TCR, significant_TCR_index, test_index )
        Auc_i = constrained_pred(HLA, HLA_index, CMV, test_index, pred_probs)
        np.append(temp_res, Auc_i)
        
    res_dict[HLA_index] = temp_res
    return


if __name__ == "__main__":
    data_path = "/home/hyo/mystuff/hutch_project/data"
    res_path =  "/home/hyo/mystuff/hutch_project/data/result"
    res_file_name = 'split_all_res.txt'

    type1_error = 0.001
    TCR, HLA, CMV, split_mat, mid_HLA_index, train_index_mat, test_index_mat = read_in_data(data_path)

    manager = multiprocessing.Manager()
    res_dict =manager.dict()
    pss = []

    for pid in range(3):
        p = multiprocessing.Process(target=work_wrapper, args=[TCR, HLA, CMV, train_index_mat[pid],
                                                               test_index_mat[pid], type1_error, pid, split_mat, res_dict])
        p.start()
        pss.append(p)
    for proc in pss:
        proc.join()

    for key in res_dict.keys():
        with open(os.path.join(res_path, res_file_name), 'a') as file1:
            file1.writelines( str(key)+","+str(res_dict[key])+"\n") 





