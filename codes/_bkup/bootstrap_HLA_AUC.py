import multiprocessing
from scipy import io
import pandas as pd
import numpy as np
import os
import time
import pickle
from scipy.stats import fisher_exact
from sklearn import metrics


def filter_p_vals(cutoff, p_vals):
    p_val_res = np.zeros(len(p_vals))
    index_res = np.zeros(len(p_vals), int)
    res_index = 0
    for i in range(len(p_vals)):
        if p_vals[i] <= cutoff:
            p_val_res[res_index] = p_vals[i]
            index_res[res_index] = int(i)
            res_index += 1
    
    return(index_res[0:res_index])


def make_pred(TCR, important_index, test_individuals):
    #HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
    #training_TCR = TCR[:, HLA_sharing_individuals]
    testing_TCR = TCR[:, test_individuals]
    important_TCRs = [ testing_TCR[i] for i in important_index ]
    TCR_of_interest = np.sum(important_TCRs, axis=0, dtype=int) # x_i for individual i
    Total_TCR = np.sum(testing_TCR, axis=0, dtype=int) # d_i
    prediction_probs = TCR_of_interest/Total_TCR # x_i / d_i
    return(prediction_probs)



def get_fisher_data(training_index, TCR):
    cardi = len(TCR)
    oneone_vec = np.zeros(cardi)
    onezero_vec = np.zeros(cardi)
    zeroone_vec = np.zeros(cardi)
    zerozero_vec = np.zeros(cardi)
    res_index = 0
    
    if training_index.size == 0:
        return np.array([oneone_vec, onezero_vec, zeroone_vec, zerozero_vec])
    
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
    return np.array([oneone_vec, onezero_vec, zeroone_vec, zerozero_vec])


def read_in_data(data_path):

    TCR_sparsematrix = io.mmread( os.path.join(data_path, "TCR.txt"))
    HLA_sparsematrix = io.mmread( os.path.join(data_path, "HLA.txt") )
    df = pd.read_csv(os.path.join(data_path, "CMV.txt") )
    TCR = TCR_sparsematrix.toarray()
    HLA = HLA_sparsematrix.toarray()
    CMV = df['x']
    CMV = CMV.to_numpy()
    return TCR, HLA, CMV



def create_train_condition_HLA(TCR, HLA, CMV, HLA_index):
    HLA_sharing_individuals = np.where(HLA[HLA_index,]  > 0)[0]
    train_len = int(len(HLA_sharing_individuals)*0.8)
    train_indexes =  np.random.choice(HLA_sharing_individuals, train_len, False)
    test_indexes = [element for element in HLA_sharing_individuals if element not in train_indexes]
    return train_indexes, test_indexes


def asso_test(mat_11, mat_10, mat_01, mat_00):
    cardi = 1098738
    res = np.zeros(cardi)
    for i in range(cardi):
        asso_i = fisher_exact( [ [ mat_11[i], mat_10[i] ], [ mat_01[i], mat_00[i] ]], alternative="greater")
        res[i] = asso_i.pvalue
    #output_dict[HLA_index] = res
    return res


def setting_targets(input_array_indexes, desired_cores):
    if input_array_indexes.size == 0:
        targets = np.array([])
        return targets
    
    if 0 < input_array_indexes.size <= desired_cores:
        num_workers = input_array_indexes.size
        targets = np.array_split(input_array_indexes, num_workers)
        return targets

    if input_array_indexes.size > desired_cores:
        num_workers = desired_cores
        targets = np.array_split(input_array_indexes, num_workers)
        return targets



def work_wrapper(TCR, HLA, CMV, targets, type1_error, num_bootstrap_samples, res_dict):
    for target_HLA in targets:
        AUCs = []
        previous_indexes = np.array([])
        cardi = len(TCR)
        oneone_vec = np.zeros(cardi)
        onezero_vec = np.zeros(cardi)
        zeroone_vec = np.zeros(cardi)
        zerozero_vec = np.zeros(cardi)
        
        
        for sample_i in range(num_bootstrap_samples):
            train_indexes, test_indexes = create_train_condition_HLA(TCR, HLA, CMV, target_HLA)
            indexes_to_remove = np.setdiff1d(previous_indexes, train_indexes)
            indexes_to_add = np.setdiff1d(train_indexes, previous_indexes)

            del_vec = get_fisher_data(indexes_to_remove, TCR)
            add_vec = get_fisher_data(indexes_to_add, TCR)

            oneone_vec = oneone_vec - del_vec[0] + add_vec[0]
            onezero_vec = onezero_vec - del_vec[1] + add_vec[1]
            zeroone_vec = zeroone_vec - del_vec[2] + add_vec[2]
            zerozero_vec = zerozero_vec - del_vec[3] + add_vec[3]

            p_vals = asso_test(oneone_vec, onezero_vec,zeroone_vec, zerozero_vec)
            significant_TCR_index = filter_p_vals(type1_error, p_vals)
            pred_probs = make_pred(TCR, significant_TCR_index, test_indexes )
            Auc_i = metrics.roc_auc_score( CMV[ test_indexes ], pred_probs)
            AUCs.append(Auc_i)
            previous_indexes = train_indexes
            
            
        res_dict[target_HLA] = AUCs
    return



if __name__ == "__main__":
    data_path = "/home/hyo/mystuff/hutch_project/data"
    bootstrap_n = 10
    type1_error = 0.001
    HLA_shared_by = 60
    desired_cores = 12 
    
    
    TCR,HLA,CMV = read_in_data(data_path)
    target_HLAs = np.where( HLA.sum(axis = 1) > HLA_shared_by )[0]
    if len(target_HLAs) <= desired_cores:
        num_workers = len(target_HLAs)
    else:
        num_workers = desired_cores

    work_targets  = np.array_split(target_HLAs, num_workers)
    manager = multiprocessing.Manager()
    res_dict =manager.dict()

    pss = []
    start = time.time()
    for pid in range(num_workers):
        p = multiprocessing.Process(target=work_wrapper, args=[TCR, HLA, 
                                                               CMV, work_targets[pid],
                                                               type1_error, bootstrap_n, res_dict])
        p.start()
        pss.append(p)
    for proc in pss:
        proc.join()
    end = time.time()
    print(end-start)
    mydict = dict()
    for key in res_dict.keys():
        mydict[key] = res_dict[key]
    with open('Feb6_AUC.pkl', 'wb') as f:
        pickle.dump(mydict, f)
            

    
