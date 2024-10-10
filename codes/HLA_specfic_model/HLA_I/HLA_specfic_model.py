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
    AA_index = np.array(pd.read_csv(os.path.join(data_path, "AA_HLA_index.csv"),header=None,dtype=int )[0])
    return TCR, HLA, CMV, train_index, test_index, AA_index



def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    return res


def create_train_condition_HLA(TCR, HLA, CMV, HLA_index):
    # notice we are indexing 1 here since we are testing condition for a 
    # 2-dimesional array, it returns two arrays
    HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[1]
    train_len = int(len(HLA_sharing_individuals)*0.8)
    train_indexes =  np.random.choice(HLA_sharing_individuals, train_len, False)
    test_indexes = [element for element in HLA_sharing_individuals if element not in train_indexes]
    return train_indexes, test_indexes


def constrained_pred(HLA_i, CMV, test_indexes, pred_probs):
    position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
    HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
    position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
    try:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
    except ValueError:
        Auc_i = -1
    return Auc_i
    

def work_wrapper(TCR, HLA, CMV, AA_index, train_index, test_indexes, type1_error, res_dict, pid, res_path):
    cardi = len(TCR)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    target_HLA = AA_index[pid]
    HLA_sharing_individuals = np.where(HLA[target_HLA,]  > 0)[0]
    train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
    oneone_vec, onezero_vec, zeroone_vec, zerozero_vec = get_fisher_data(train_indexes, TCR)
    p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
    significant_TCR_index, signi_TCR_pvals = filter_p_vals(type1_error, p_vals)
    pred_probs = make_pred(TCR, significant_TCR_index, test_indexes )
    Auc_i = constrained_pred(target_HLA, CMV, test_indexes, pred_probs)
    res_dict[target_HLA] = Auc_i
    return


# Parallizing the HLA.
# The prediction is constrained to the test individuals that have this particular HLA.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The group_id for HLA_I allele indexes')
    parser.add_argument('--group_id', type=int, help='group_id in integer format')
    input_args = parser.parse_args()
    nums = input_args.sample_id
    start = int(nums*10)
    if start == 80:
        end = 85
    else:
        end = start + 10

    os.chdir("./")
    data_path ="../../intermediate_files/"
    res_path ="./"
    res_file_name = '{}_group_res.txt'.format(start)


    type1_error = 0.001
    TCR, HLA, CMV, train_index, test_index, AA_index = read_in_data(data_path)
    manager = multiprocessing.Manager()
    res_dict =manager.dict()
    pss = []

    for pid in range(start, end):
        p = multiprocessing.Process(target=work_wrapper, args=[TCR, HLA, CMV, AA_index,
                                                               train_index, test_index, type1_error, res_dict, pid, res_path])
        p.start()
        pss.append(p)
    for proc in pss:
        proc.join()
    
    for key in res_dict.keys():
        with open(os.path.join(res_path, res_file_name), 'a') as file1:
            file1.writelines( str(key)+","+str(res_dict[key])+"\n") 

