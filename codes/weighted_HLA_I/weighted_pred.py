#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import os
import argparse
import multiprocessing
import numpy as np
from scipy import io
import pandas as pd
from sklearn import metrics
from scipy.stats import chi2
from datetime import datetime
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning
import warnings
warnings.simplefilter('ignore', PerfectSeparationWarning)
warnings.filterwarnings('ignore') 


def drop_df_nas(input_df, flag):
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
    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    train_index = np.array(pd.read_csv(os.path.join(data_path, "train_index.csv"),header=None,dtype=int )[0])
    test_index = np.array(pd.read_csv(os.path.join(data_path, "test_index.csv"),header=None,dtype=int )[0])
    HLA_index_mat = pd.read_csv( os.path.join(data_path, "HLA_I_neighbor_index_mat.csv"))
    HLA_index_mat = drop_df_nas(HLA_index_mat, 1)
    neigbor_weights_mat = pd.read_csv( os.path.join(data_path, "HLA_I_weights_mat.csv"))
    neigbor_weights_mat = drop_df_nas(neigbor_weights_mat, 0)
    AA_index = np.array(pd.read_csv(os.path.join(data_path, "HLA_I_index.csv"),header=None,dtype=int )[0])
    return TCR, HLA, CMV, train_index, test_index, HLA_index_mat, neigbor_weights_mat, AA_index


def get_p_vals(TCR_i, CMV, log_rd, train_index, train_weights):
    X = np.column_stack((CMV[train_index], log_rd[train_index].reshape(-1,1), TCR_i[train_index] ))
    X = pd.DataFrame(X)
    X.columns =["CMV", "log_rd", "TCR_i"]
    if ((len(set(X["CMV"]))==1) or (len(set(X["TCR_i"]))==1)):
        lr_pvalue = np.nan
    else:
        results = smf.glm("TCR_i~log_rd+CMV", family=sm.families.Binomial(), data=X, freq_weights = train_weights).fit()
        if results.params[2] < 0:
            lr_pvalue = 1
            return lr_pvalue

        results_res = smf.glm("TCR_i~log_rd", family=sm.families.Binomial(), data=X, freq_weights = train_weights).fit()
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

    HLA_train_slice_2 = HLA_slice_1[:, train_indexes]
    weight_mat = HLA_train_slice_2.transpose() * weights_values
    train_weights = np.max(weight_mat, axis=1)
    return train_indexes, train_weights


def constrained_pred(target_HLA, CMV, test_indexes, pred_probs):
    test_indexes = np.array(test_indexes)
    AUC = []
    for HLA_i in target_HLA:
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


def work_for_HLA_i(HLA_index_mat, neigbor_weights_mat, HLA, TCR, HLA_index, CMV, log_rd, train_index, test_index, 
                   HLA_I_index, pval_path, auc_path):
    target_HLA = HLA_I_index[HLA_index]
    train_indexes, train_weights = train_test_split(HLA, HLA_index_mat[HLA_index], train_index, neigbor_weights_mat[HLA_index])

    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%H:%M:%S")
    print(current_time)
    print("start computing p-values")

    important_indexes = []
    TCR_pvals = []
    for i in range(len(TCR)):
        p_val_i = get_p_vals(TCR[i], CMV, log_rd, train_indexes, train_weights)
        TCR_pvals.append(p_val_i)
        if p_val_i <= 0.001:
            important_indexes.append(i)
            
    important_indexes = np.array(important_indexes)

    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%H:%M:%S")
    print(current_time)
    print("finished computing p-values")

    pval_filename = os.path.join(pval_path, f"pvalues_hla_index_{target_HLA}.csv")
    df_pvals = pd.DataFrame(TCR_pvals, 
                            columns=['pval'])
    df_pvals.to_csv(pval_filename, index=False, na_rep="NA")
    
    if important_indexes.size == 0:
        pred_probs = np.zeros(len(test_index))
    else:
        pred_probs = TCR[important_indexes][:,test_index].sum(axis=0) / TCR[:, test_index].sum(axis=0)

    Aucs = constrained_pred(HLA_index_mat[HLA_index], CMV, test_index, pred_probs)

    auc_filename = os.path.join(auc_path, f"auc_hla_index_{target_HLA}.csv")   
    df_auc = pd.DataFrame(Aucs, 
                          columns=['AUC'])
    df_auc.to_csv(auc_filename, index=False, na_rep="NA")

    return


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='prediction based on weighted K nearest neighbors')
    parser.add_argument('--hla_id', type=int, help='the relative index of HLA (among HLA-Is)')
    input_args = parser.parse_args()
    hla_id = input_args.hla_id

    print(hla_id)  

    pval_path = "./results/pvalues"
    auc_path = "./results/res_aucs"
    os.makedirs(pval_path, exist_ok=True)
    os.makedirs(auc_path, exist_ok=True)

    data_path = "../intermediate_files"
    TCR, HLA, CMV, train_index, test_index, HLA_index_mat, neigbor_weights_mat, HLA_I_index = read_in_data(data_path)
    log_rd = np.log(np.sum(TCR,0))

    work_for_HLA_i(HLA_index_mat, neigbor_weights_mat, HLA, TCR, hla_id, CMV, log_rd,
                   train_index, test_index, HLA_I_index, pval_path, auc_path)
