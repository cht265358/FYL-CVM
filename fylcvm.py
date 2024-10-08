#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import matplotlib.pyplot as plt
from scipy.optimize import minimize
#svm stands for site_variable_matrix

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
        #print(self.basicclusterenergy)
        print("finish mapping")
    
    def compute_site_variable_matrix(self,site_variable_input,composition,T):         #extend it to larger system if expression is known
        svm=np.ones((self.component,self.clustersize))
        for i in range(self.clustersize-1):
            svm[1][i]=np.exp(site_variable_input[i])                                  #the input is the log of site variable
        nominator=0
        denominator=0
        shape1=[self.component]*(self.clustersize-1)
        for i,j,k in itertools.product(range(self.component), repeat=3):
            nominator+=(composition[1]-(1/self.clustersize)*(i+j+k+0))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][0]/T)
            denominator+=(-composition[1]+(1/self.clustersize)*(i+j+k+1))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][1]/T)
        svm[1][3]=nominator/denominator
        #print(svm[1][3])
        return svm
        
    def compute_partition_function(self,site_variable_input,composition,T):
        svm=self.compute_site_variable_matrix(site_variable_input,composition,T)      #use A as reference mu(A)=0
        if self.clustersize==4:
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp(-self.basicclusterenergy[i][j][k][l]/T)
            return partition
        else:
            print("unfinished")
            return 1;

    def compute_basic_cluster_prob(self,partition,site_variable_input,composition,T):
        svm=self.compute_site_variable_matrix(site_variable_input,composition,T)             #svm is site potential matrix
        basic_cluster_prob=np.zeros(self.shape)
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            basic_cluster_prob[i][j][k][l]=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp(-self.basicclusterenergy[i][j][k][l]/T)/partition
        return basic_cluster_prob

    def compute_total_energy(self,site_variable_input,T,component_comp):    #site_variable_input is 1d for binary for now
        svm=self.compute_site_variable_matrix(site_variable_input,component_comp,T) 
        partition=self.compute_partition_function(site_variable_input,component_comp,T)
        prob=self.compute_basic_cluster_prob(partition,site_variable_input,component_comp,T)
        #print(prob)
        #print(partition)
        if partition<0:
            return 1000.0
        #two body cluster
        #print("two body cluster")
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
        #point cluster  
        #print("point cluster")          
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=T*np.log(partition)
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        return total_energy
                
    def optimize_free_energy(self,T,component_comp):         #main optimization function, optimize site variables at given T and composition
        #initial_guess=np.ones((self.clustersize-1)*(self.component-1))
        initial_guess2=np.array([[5,5,5],[5,5,-5],[2,2,0],[10,10,10]])
        guess=np.array([10,10,10])
        options={
            'initial_simplex':initial_guess2
        }
        positive=((-20,20),(-20,20),(-20,20))
        #result=minimize(CVM.compute_total_energy,guess,method='Nelder-Mead',args=(T,component_comp),bounds=positive,tol=1.0e-6)
        result=minimize(CVM.compute_total_energy,guess,method='Nelder-Mead',options=options,args=(T,component_comp),bounds=positive,tol=1.0e-6)
        return result
    
    def output_optimization(self,result,T,component_comp):
        svm=self.compute_site_variable_matrix(result.x,component_comp,T) 
        print("free energy is "+str(result.fun)+" site variable is "+str(svm[1])+" at temperature "+str(T))

    def print_output(self):
        print("hello world")

if __name__ == '__main__':
    print("FYLCVM code")
    E=[0,-3,-4,-3,0]
    CVM=FYLCVM(2,4)
    CVM.map_basic_cluster_energy(E)
    composition=np.array([0.4,0.6])
    #result=CVM.optimize_free_energy(1.0,composition)
    
    '''
    Tspace=np.linspace(1.0,2.5,61)
    for T in Tspace:
        svm=CVM.compute_site_variable_matrix(np.array([np.log(6.95524423),np.log(6.95524423),np.log(6.95524423)]),composition,T)
        print("T is "+str(T)+"  "+str(svm[1][3]))
    
    #test first with fixed mu
    '''
    count=0
    Tspace=np.linspace(1.0,2.0,41)
    Fspace=np.zeros(41)
    for T in Tspace:
        result=CVM.optimize_free_energy(T,composition)
        CVM.output_optimization(result,T,composition)
        Fspace[count]=result["fun"]
        count+=1
    plt.plot(Tspace,Fspace)
    plt.show()
    #print("free energy is "+str(F)+"entropy is "+str(S)+" at T="+str(T))
    