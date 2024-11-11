#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit

#this code can do the following tasks
#maps the energies of basic clusters
#propose structures to compute basic cluster energies

def TG(T,a,b,c,d,e):
    return a+b*T+c*T**2+d*T*np.log(T)+e/T

def TG_easy(T,params):
    return params[0]+params[1]*T+params[2]*T**2+params[3]*T*np.log(T)+params[4]/T

class mapper:
    def __init__(self,lattice,filename):
        self.TF=np.loadtxt(filename)
        self.lattice=lattice
    
    def fit(self):
        print("place holder for overwritting")

class F_FCC(mapper):
    def __init__(self,lattice,filename="energy.in"):
        super().__init__(self,lattice,filename)
        self.Tmat=self.TF[:,0]
        self.mAAAA=self.TF[:,1]
        self.mBBBB=self.TF[:,5]
        self.mAAAB=self.TF[:,2]-0.75*self.mAAAA-0.25*self.mBBBB
        self.mAABB=self.TF[:,3]-0.5*self.mAAAA-0.5*self.mBBBB
        self.mABBB=self.TF[:,4]-0.25*self.mAAAA-0.75*self.mBBBB

    def AAAB(self,T):
        fAAAB, covariance = curve_fit(TG,self.Tmat,self.mAAAB)
        return TG_easy(T,fAAAB)
    
    def AABB(self,T):
        fAABB, covariance = curve_fit(TG,self.Tmat,self.mAABB)
        return TG_easy(T,fAABB)
    
    def ABBB(self,T):
        fABBB, covariance = curve_fit(TG,self.Tmat,self.mABBB)
        return TG_easy(T,fABBB)


if __name__ == '__main__':
    print("This is mapper")