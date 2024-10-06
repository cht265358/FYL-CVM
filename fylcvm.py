#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
#from scipy.optimize import basinhopping    algorithm placeholder
#spm stands for site_potential_matrix

class FYLCVM:       #base class for FCC
    #
    def __init__(self,component,maxsize):                #initialize all variables
        print("hello world")
        kB=1.380649e-23  #boltzman constant, unit:J/K
        N=6.02e23;       #number of particles, the size of the system, here we use 1 mole
        self.lattice="FCC"
        self.component=component
        self.clustersize=maxsize
        self.shape=[self.component]*self.clustersize    #ignore using basic cluster with different size first

    def map_basic_cluster_energy(self,E):     #map the basic cluster energy to high dimensional array
        self.basicclusterenergy=np.zeros(self.shape)                  #use high dimensional array, input by hand now
        if self.clustersize==4:
            for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
                self.basicclusterenergy[i][j][k][l]=E[i+k+j+l]
        print(self.basicclusterenergy)
        print("finish mapping")
    
    def compute_partition_function(self,site_potential):
        spm=np.append(np.zeros(self.clustersize),site_potential)      #use A as reference mu(A)=0
        if self.clustersize==4:
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=np.exp(spm[i][0]+spm[j][1]+spm[k][2]+spm[l][3]-self.basicclusterenergy[i][j][k][l])
            return partition
        else:
            print("unfinished")
            return 1;

    def compute_basic_cluster_prob(self,partition,site_potential):
        spm=np.append(np.zeros(self.clustersize),site_potential)               #spm is site potential matrix
        basic_cluster_prob=np.zeros(self.shape)
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            basic_cluster_prob[i][j][k][l]=np.exp(spm[i][0]+spm[j][1]+spm[k][2]+spm[l][3]-self.basicclusterenergy[i][j][k][l])/partition

    def compute_total_energy(self,site_potential,prob,T):
        partition=self.compute_partition_function(site_potential)
        #two body cluster
        two_body_energy=0
        for i in range(self.clustersize):                 #4*3*2*2 two body cluster
            for j in range(i+1,self.clustersize):
                for k,l in itertools.product(range(self.component), repeat=2):     #decorate this line later
                    position_type_matrix=np.zeros((2,2))
                    position_type_matrix[0][0]=i
                    position_type_matrix[0][1]=j
                    position_type_matrix[1][0]=k
                    position_type_matrix[1][1]=l
                    probij=utility.get_cluster_prob(prob,position_type_matrix)
                    two_body_energy+=probij*np.log(probij)*T
        pointenergy=0
        point_prob=np.zeros(self.maxsize)
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                probi=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=probi*np.log(probi)*T
        totalpotential=0
        for 

                
        

    def optimize_free_energy(self):         #main optimization function
        self.optimize()         #non-linear optimization utility function

    def optimize(self):
        print("hello world")

if __name__ == '__main__':
    print("FYLCVM code")
    E=[0,-3,-4,-3,0]
    CVM=FYLCVM(2,4)
    CVM.map_basic_cluster_energy(E)
    mu=np.array([1,1,-1,-1])
    partition=CVM.compute_partition_function(mu)
    print(partition)