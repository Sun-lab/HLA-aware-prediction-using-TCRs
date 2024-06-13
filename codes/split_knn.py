#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 19:34:50 2024

@author: hyo
"""

import os
os.environ["OMP_NUM_THREADS"] = "16" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "16" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "16" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "16" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "16" # export NUMEXPR_NUM_THREADS=6
import argparse
import multiprocessing
import numpy as np
from scipy import io
import pandas as pd
from sklearn import metrics
from scipy.stats import chi2
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning
import warnings
warnings.simplefilter('ignore', PerfectSeparationWarning)
warnings.filterwarnings('ignore') 


def drop_df_nas(input_df, flag):
    if flag == 1:
        input_df = input_df - 1
    input_df = input_df.to_numpy()
    test = input_df.tolist()
    res = []
    for sub_array in test:
        if flag == 1:
            res.append([int(x) for x in sub_array if str(x) != 'nan'])
        else:
            res.append([float(x) for x in sub_array if str(x) != 'nan'])
    return res
    

def read_in_data(data_path):
    exception_index = [88]
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    HLA_index_mat = pd.read_csv( os.path.join(data_path, "index_mat.csv"))
    HLA_index_mat = drop_df_nas(HLA_index_mat, 1)
    neigbor_weights_mat = pd.read_csv( os.path.join(data_path, "weights_mat.csv"))
    neigbor_weights_mat = drop_df_nas(neigbor_weights_mat, 0)
    mid_HLA_index = np.array(pd.read_csv(os.path.join(data_path, "mid_HLA_index.csv"), header=None, dtype=int)[0])
    train_index_mat = np.array(pd.read_csv(os.path.join(data_path, "split_train_index.csv"), header=None, dtype=int))
    test_index_mat = np.array(pd.read_csv(os.path.join(data_path, "split_test_index.csv"), header=None, dtype=int))
    split_mat = np.array(pd.read_csv(os.path.join(data_path, "split_mat.csv"), header=None, dtype=int))
    new_HLA = np.vstack((HLA, HLA[exception_index]))


    return TCR, new_HLA, CMV, train_index_mat, test_index_mat, HLA_index_mat, neigbor_weights_mat, mid_HLA_index, split_mat


def get_p_vals(TCR_i, CMV, log_rd, train_index, train_weights):
    X = np.column_stack((TCR_i[train_index], log_rd[train_index].reshape(-1,1), CMV[train_index] ))
    X = pd.DataFrame(X)
    X.columns =["Y", "log_rd", "CMV"]
    results = smf.glm("Y~log_rd+CMV", family=sm.families.Binomial(), data=X, freq_weights = train_weights).fit()
    results_res = smf.glm("Y~log_rd", family=sm.families.Binomial(), data=X, freq_weights = train_weights).fit()
    llf_full = results.llf
    llf_restr = results_res.llf
    df_full = results.df_resid 
    df_restr = results_res.df_resid 
    lrdf = (df_restr - df_full)
    lrstat = -2*(llf_restr - llf_full)
    lr_pvalue = chi2.sf(lrstat, df=lrdf)
    return lr_pvalue


def train_test_split(HLA, neigbhor_index, train_index, weights_values):
    HLA_slice_1 = HLA[neigbhor_index, ]
    colsums_HLA_slice = HLA_slice_1.sum(axis=0)
    whole_individuals = np.where(colsums_HLA_slice > 0)[0]
    
    train_indexes = np.intersect1d(whole_individuals, train_index)
    #train_len = int(len(whole_individuals)*0.8)
    #train_index =  np.random.choice(whole_individuals, train_len, False)
    #test_index = np.setdiff1d(whole_individuals, train_index)
    HLA_train_slice_2 = HLA_slice_1[:, train_indexes]
    weight_mat = HLA_train_slice_2.transpose() * weights_values
    train_weights = np.max(weight_mat, axis=1)
    return train_indexes, train_weights


def constrained_pred(HLA, target_HLA, CMV, test_indexes, pred_probs):
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
    return AUC


def work_for_HLA_i(HLA_index_mat, neigbor_weights_mat, HLA, TCR, HLA_index, CMV, log_rd, res_dict, train_index, test_index, HLA_I_index ):
    target_HLA = HLA_I_index[HLA_index]
    train_indexes, train_weights = train_test_split(HLA, HLA_index_mat[HLA_index], train_index, neigbor_weights_mat[HLA_index])
    important_indexes = []
    TCR_pvals = []
    for i in range(len(TCR)):
        p_val_i = get_p_vals(TCR[i], CMV, log_rd, train_indexes, train_weights)
        TCR_pvals.append(p_val_i)
        if p_val_i <= 0.001:
            important_indexes.append(i)
            
    pval_array = np.array(TCR_pvals)
    important_indexes = np.array(important_indexes)
    res_name1 = os.path.join("/home/fding/weighted_pred_HLAII/results", "{}TCR_pvals.txt".format(target_HLA))
    np.savetxt(res_name1, pval_array)
    
    if important_indexes.size == 0:
        pred_probs = np.zeros(len(test_index))
    else:
        pred_probs = TCR[important_indexes][:,test_index].sum(axis=0) / TCR[:, test_index].sum(axis=0)
        
    Aucs = constrained_pred(HLA_index_mat[HLA_index], CMV, test_index, pred_probs)
    res_dict[HLA_index] = {"cluster" : HLA_index_mat[HLA_index], "AUC" : Aucs }

    return


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='sample_id plus 6, repeat 10 times')
    parser.add_argument('--sample_id', type=int, help='sample_id in integer format')
    input_args = parser.parse_args()
    nums = input_args.sample_id
    start = int(nums*10)
    if start == 130:
        end = 135
    else:
        end = start + 10
        
        
    res_path = "/home/fding/weighted_pred_HLAII/results"
    res_file_name = '{}_group_res.txt'.format(start)

    data_path = "/home/fding/weighted_pred_HLAII/data_files"
    TCR, new_HLA, CMV, train_index_mat, test_index_mat, HLA_index_mat, neigbor_weights_mat, mid_HLA_index, split_mat = read_in_data(data_path)
    log_rd = np.log(np.sum(TCR,0))
    manager = multiprocessing.Manager()
    res_dict =manager.dict()
    pss = []

    for pid in range(start, end):
        p = multiprocessing.Process(target=work_for_HLA_i, args=[HLA_index_mat, neigbor_weights_mat, HLA, TCR, pid, CMV, log_rd,
                                                                 res_dict, train_index, test_index, HLA_II_index])
        p.start()
        pss.append(p)
    for proc in pss:
        proc.join()
    
    for key in res_dict.keys():
        with open(os.path.join(res_path, res_file_name), 'a') as file1:
            file1.writelines(str(key)+","+str(res_dict[key])+"\n") 
    

