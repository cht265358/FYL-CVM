#!/usr/bin/python

import numpy as np
import itertools

#This is the utility function file include all the "dirty" funtions

def get_cluster_prob(high_d_matrix,position_type_matrix):     #get cluster probablity from high dimensional matrix
    if np.sum(high_d_matrix)==1.0:
        shape=np.shape(high_d_matrix)
        prob=0
        if np.shape(position_type_matrix)[1]==1:     #compute point cluster probability
            location=position_type_matrix[0][0]
            atom_type=position_type_matrix[1][0]
            for idx in itertools.product(*[range(s) for s in shape]):
                if idx[location]==atom_type:
                    prob+=high_d_matrix[idx]
        if np.shape(position_type_matrix)[1]==2:     #compute binary cluster probability
            location=position_type_matrix[0]
            atom_type=position_type_matrix[1]
            for idx in itertools.product(*[range(s) for s in shape]):
                if idx[location[0]]==atom_type[0] and idx[location[1]]==atom_type[1]:
                    prob+=high_d_matrix[idx]
        return prob

    else:
        print("error, total probability is not equal to 1")
        return 0
    
if __name__ == '__main__':
    arr = np.ones([2,2,2,2])
    arr = arr/np.sum(arr)
    print(np.sum(arr))
    position_type_matrix=np.array([[2,3],[1,1]],ndmin=2)
    print(get_cluster_prob(arr,position_type_matrix))