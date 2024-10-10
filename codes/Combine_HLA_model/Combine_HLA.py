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
    return TCR, new_HLA, CMV, train_index, test_index





def predict_test_individual(test_index, potential_HLAs, HLA, TCR, HLA_specifc_TCRs, base_index ):
    pred_probs = []
    
    for subj_index, ind_index in enumerate(test_index):
        carry_flag = HLA[:,ind_index][potential_HLAs]
        if carry_flag == 1:
            individual_TCR_index = np.concatenate([base_index, HLA_specifc_TCRs])
            individual_TCR_index = np.unique(individual_TCR_index)    
            
        if carry_flag == 0:
            individual_TCR_index = base_index
        
        prop = make_pred(TCR, individual_TCR_index, ind_index )
        pred_probs.append(prop)
        
    return np.array(pred_probs)


def constrained_pred(HLA, HLA_i, CMV, test_indexes, pred_probs):
    test_indexes = np.array(test_indexes)
    position_index = np.arange( len(HLA[HLA_i][test_indexes]) )
    HLA_i_test_index = test_indexes[ HLA[HLA_i][test_indexes] != 0 ]
    position_index_pred = position_index[ HLA[HLA_i][test_indexes] != 0  ]
    try:
        Auc_i = metrics.roc_auc_score( CMV[ HLA_i_test_index ], pred_probs[position_index_pred] )
    except ValueError:
        Auc_i = -1
    return Auc_i



if __name__ == "__main__":
    TCR, HLA, CMV, train_index, test_index = read_in_data("/home/fding/specific_HLA_II_pred/data_files")
    all_tcr_pvals = np.loadtxt("/home/fding/all_sample_pred/results/allTCR_pvals.txt")
    base_tcr_pool,_ = filter_p_vals(0.001, all_tcr_pvals)
    HLA_test_index = np.arange(220)
    AUCS3 = []
    AUCS4 = []
    AUCS5 = []
    # this for loop reads in all the HLA pvals and create the HLA index <-> significant TCR index map
    for i in HLA_test_index:
        pvals = np.loadtxt("/home/fding/pvals/three/{}TCR_pvals.txt".format(i))

        important_index3, pvalss3 = filter_p_vals(0.001, pvals)
        important_index4, pvalss4 = filter_p_vals(0.0001, pvals)
        important_index5, pvalss5 = filter_p_vals(0.00001, pvals)

    
        pred_probs3 = predict_test_individual(test_index, i, HLA, TCR, important_index3, base_tcr_pool )
        pred_probs4 = predict_test_individual(test_index, i, HLA, TCR, important_index4, base_tcr_pool )
        pred_probs5 = predict_test_individual(test_index, i, HLA, TCR, important_index5, base_tcr_pool )

    

        AUCS3.append(constrained_pred(HLA, i, CMV, test_index, pred_probs3))
        AUCS4.append(constrained_pred(HLA, i, CMV, test_index, pred_probs4))
        AUCS5.append(constrained_pred(HLA, i, CMV, test_index, pred_probs5))
    
    df345 = pd.DataFrame({'AUC3':AUCS3, 'AUC4':AUCS4, 'AUC5':AUCS5})
    df345.to_csv("diff_aucs.csv", sep='\t', encoding='utf-8', index=False, header=False)
    

