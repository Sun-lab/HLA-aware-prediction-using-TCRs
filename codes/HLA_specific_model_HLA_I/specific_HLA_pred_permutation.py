#!/usr/bin/env python3
# -*- coding: utf-8 -*-




import multiprocessing
import numpy as np
from scipy import io
import pandas as pd
import os
from scipy.stats import fisher_exact
from sklearn import metrics
import random
from datetime import datetime
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
    #HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
    #training_TCR = TCR[:, HLA_sharing_individuals]
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
    

def work_wrapper(TCR, HLA, CMV, AA_index, train_index, test_indexes, type1_error, pid):
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
#    res_name = os.path.join(res_path, "{}TCR_index.txt".format(target_HLA))
#    res_name1 = os.path.join(res_path, "{}TCR_pvals.txt".format(target_HLA))

#    np.savetxt(res_name, significant_TCR_index)
#    np.savetxt(res_name1, p_vals)

    pred_probs = make_pred(TCR, significant_TCR_index, test_indexes )
    Auc_i = constrained_pred(target_HLA, CMV, test_indexes, pred_probs)
#    res_dict[target_HLA] = Auc_i
    return Auc_i

'''
AA_index[0]

train_size = []
test_size = []


for target_HLA in AA_index:
    HLA_sharing_individuals = np.where(HLA[target_HLA,]  > 0)[0]
    train_indexes = np.intersect1d(train_index, HLA_sharing_individuals)
    train_size.append(len(train_indexes))
    test_size.append(len(test_index[ HLA[target_HLA][test_index] != 0 ]))

train_size = np.array(train_size)
test_size = np.array(test_size)

df = pd.DataFrame({'HLA_index': AA_index+1, 'train_size_specfic': train_size, 'test_size_specfic': test_size})
df.to_csv('/home/hyo/mystuff/hutch_project/specific_size.csv', index=False)
#AA_index
'''

# no clustering, for each HLA, we build the model using all the individuals that carry this HLA.
# Recommanding parallize the HLA.
# The prediction is constrained to the test individuals that have this particular HLA.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='specific HLA prediction permutation p-values by chunk')
    parser.add_argument('--hla_id', type=int, help='the index of HLA (among HLA-Is) to run permutation for in integer format')
    parser.add_argument('--chunk_id', type=int, help='the id of the chunk of permutation job submission')
    input_args = parser.parse_args()
    hla_id = input_args.hla_id
    chunk_id = input_args.chunk_id

    print(hla_id)
    print(chunk_id)

    thickness = 100

    data_path ="./data_files"
    res_path ="./results/permutation_aucs/hla_index_"+str(hla_id)
    os.makedirs(res_path, exist_ok=True)

    res_file_name = str(hla_id)+"_"+str(chunk_id)+"_chunk_res.csv"

    # load corresponding random seeds
    df_seeds = pd.read_csv(os.path.join(data_path, "random_seed_100.txt"), 
                           sep=" ", header=None)
    seeds_100 = df_seeds[0].tolist()
    cur_seed = seeds_100[chunk_id]

    TCR,HLA_ori,CMV,train_index,test_index,AA_index = read_in_data(data_path)

    random.seed(cur_seed)

    type1_error = 0.001

    res_list = []

    # pid = 23

    current_datetime = datetime.now()
    current_time = current_datetime.strftime("%H:%M:%S")
    print(current_time)
    print("start permutation")

    train_index_list = train_index.tolist()

    for i in range(thickness):
        permute_indexes = random.sample(train_index_list, 
                                        len(train_index_list))
        HLA = HLA_ori.copy()
        HLA[:, permute_indexes] = HLA_ori[:, train_index_list]
        # HLA_53 = HLA[53, permute_indexes].tolist()
        # HLA_ori_53 = HLA_ori[53, train_index_list].tolist()
        # sum([HLA_53[x]==HLA_ori_53[x] for x in range(len(permute_indexes))])
        # HLA_53 = HLA[53, test_index].tolist()
        # HLA_ori_53 = HLA_ori[53, test_index].tolist()
        # sum([HLA_53[x]==HLA_ori_53[x] for x in range( test_index.shape[0])])
        cur_auc = work_wrapper(TCR, HLA, CMV, AA_index, train_index, test_index, type1_error, hla_id)
        res_list += [cur_auc]
        current_datetime = datetime.now()
        current_time = current_datetime.strftime("%H:%M:%S")
        print(current_time)
        print("done with permutation "+str(i))

    df_chunk_pvalues = pd.DataFrame(res_list, 
                                    columns=['permutation_auc'])
    
    df_chunk_pvalues.to_csv(os.path.join(res_path, res_file_name), 
                            index=False)