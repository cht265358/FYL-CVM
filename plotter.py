#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.interpolate import splrep

def plot_phase_diagram(phase_boundary_list,signallist,dT,Tstart,filename="phasediagram.png",color='r'):
    print("plot phase diagram using data obtained from FYL-CVM")
    listlen=len(phase_boundary_list)
    signal_list=signallist.copy()         #protect signallist from modification
    eutecticend=[[0.4],[0.45],[589]]
    slopelist=[[0,0]]*listlen
    for i in range(listlen):                #search invariant point first
        if signal_list[i]==0:               #a phb that ends
            for j in range(i+1,listlen):
                if phase_boundary_list[i][2][-1]==phase_boundary_list[j][2][-1]:     #is invariant point
                    slopelist[i],slopelist[j]=solve_invariant(phase_boundary_list[i],phase_boundary_list[j],dT,color)
                    signal_list[i]=signal_list[j]=4
    for i in range(listlen):                #finish all endpoint
        if signal_list[i]==0:   
            slopelist[i]=solve_endpoint(phase_boundary_list[i],0.5,dT,color)
    for i in range(listlen):
        for j in range(i+1,listlen):       #search for eutectic start
            if phase_boundary_list[i][2][0]==phase_boundary_list[j][2][0] and phase_boundary_list[i][2][0]!=Tstart:
                if phase_boundary_list[i][0][0]<phase_boundary_list[j][0][0]:
                    eutecticend=solve_eutectic(phase_boundary_list[i],phase_boundary_list[j],slopelist[i],slopelist[j],dT,color)
                else:
                    #print("here is the case")
                    eutecticend=solve_eutectic(phase_boundary_list[j],phase_boundary_list[i],slopelist[j],slopelist[i],dT,color)
    for i in range(listlen):
        if signal_list[i]==2:  
            finish_eutectic(phase_boundary_list[i],eutecticend,color)
    #plt.xlim(0.1, 0.5)
    #plt.ylim(1, 2.4)
    if filename:
        plt.savefig(filename)
    #plt.show()

def solve_invariant(x_t1,x_t2,dT,color='r'):
   T=x_t1[2][-1]+dT*(x_t1[0][-1]-x_t2[0][-1])/(x_t1[0][-2]-x_t2[0][-2])
   x=0.25*(x_t1[0][-1]+x_t1[1][-1]+x_t2[0][-1]+x_t2[1][-1])
   print("invariant point estimation is x={:.3f} T={:.2f}".format(x,T))
   point=np.array([[x],[x],[T]])                #solve x based on T later, might need to create one FYLCVM class
   slope1=plot_invariant(x_t1,point,color)
   slope2=plot_invariant(x_t2,point,color)
   return slope1,slope2

def plot_invariant(xt,point,color='r'):
    xt=np.append(xt,point,axis=1)
    slope1=plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
    slope2=plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)
    return np.array([slope1,slope2])
    #plt.show()

def solve_endpoint(xt,point,dT=10,color='r'):
   Tend=xt[2][-1]+dT*((xt[0][-1]-point)/(xt[0][-2]-point))
   print("Tend is "+str(Tend))
   xt=np.append(xt,np.array([[point],[point],[Tend]]),axis=1)
   slope1=plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
   slope2=plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)
   return np.array([slope1,slope2])

def solve_eutectic(x_t1,x_t2,slope1,slope2,dT,color='r'):
    if x_t1[0][0]<x_t2[0][0]:                   #make sure xt1 on the left
        Tend=x_t1[2][0]
        x=(slope1[0]*x_t1[0][0]-slope2[1]*x_t2[1][0])/(slope1[0]-slope2[1])
        T=x_t1[2][0]+(x-x_t1[0][0])*slope1[0]
        if T<Tend-dT:
            T=Tend-dT
        x1other=x_t1[1][0]+(T-x_t1[2][0])/slope1[1]
        x2other=x_t2[0][0]+(T-x_t2[2][0])/slope2[0]
        Tspace=[Tend,T]
        xspace10=[x_t1[0][0],x]
        plt.plot(xspace10,Tspace,color)
        xspace11=[x_t1[1][0],x1other]
        plt.plot(xspace11,Tspace,color)
        xspace20=[x_t2[0][0],x2other]
        plt.plot(xspace20,Tspace,color)
        xspace21=[x_t2[1][0],x]
        plt.plot(xspace21,Tspace,color)
        xline=[x1other,x,x2other]
        Tline=[T,T,T]
        plt.plot(xline,Tline,color)
        print("T is "+str(T))
        print(np.array([[x2other],[x1other],[T]]))
        return np.array([[x2other],[x1other],[T]])

def finish_eutectic(xt,point,color='r'):
    if np.abs(xt[2][-1]-point[2][0])>1.0e-6:
        print("unequal")
        xt=np.append(xt,point,axis=1)
    else:
        #drop the last element and replace
        xt=xt[:, :-1]
        xt=np.append(xt,point,axis=1)
    slope1=plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
    slope2=plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)

def plot_phase_boundary_v1(xt,color='r'):
    sig=check_positive_v1(xt[0])
    if sig==1:
        phb=CubicSpline(xt[0],xt[1],bc_type=((2, 0.0), (1, 0.0)))
        xspace=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        plt.plot(xspace,phb(xspace),color)
        return phb(xt[0][0],1)
    elif sig==-1:
        phb=CubicSpline(1-xt[0],xt[1],bc_type=((2, 0.0), (1, 0.0)))
        xspace=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        plt.plot(xspace,phb(1-xspace),color)
        return -phb(1-xt[0][0],1)
    elif sig==0:                     #plot in yx space
        w=np.ones(len(xt[0]))
        w[-1]=10
        phb=UnivariateSpline(xt[1],xt[0],w=w)
        Tspace=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(phb(Tspace),Tspace,color)
        return 0

def check_positive_v1(x):
    length=len(x)
    if x[1]>x[0]:
        for i in range(1,length):
            if x[i]<x[i-1]:
                return 0
        return 1
    elif x[1]<x[0]:
        for i in range(1,length):
            if x[i]>x[i-1]:
                return 0
        return -1

def plot_phase_diagram_rough(phblist,filename="rough.png"):
    for phb in phblist:
        plt.plot(phb[0],phb[2])
        plt.plot(phb[1],phb[2])
    plt.savefig(filename)  
    plt.show()

def plot_scatter_rough(phblist,filename="scatter.png"):
    plt.scatter(phblist[:,0],phblist[:,1])
    plt.savefig(filename)
    plt.show()

if __name__ == '__main__':
    shell=1
    if shell:
        signallist=np.array([2,0,0,0])
        print("shell mode, input whatever you want")
        with open("e4.txt", "r") as file1:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt = eval(file1.read())
        plot_phase_diagram(xt,signallist,0.05,1,False,'g')
        print(signallist)
        with open("10.txt", "r") as file2:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt2 = eval(file2.read())
        plot_phase_diagram(xt2,signallist,0.05,1,False,'r')
        with open("e12.txt", "r") as file3:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt3 = eval(file3.read())
        plot_phase_diagram(xt3,signallist,0.05,1,False,'b')
        plt.xlim(0.1, 0.5)
        plt.ylim(1, 2.2)
        plt.xlabel("x")
        plt.ylabel("normalized T")
        plt.plot(0,0,color='r',label="($\Omega$)=0")
        plt.plot(0,0,color='g',label="($\Omega$)=4")
        plt.plot(0,0,color='b',label="($\Omega$)=12")
        plt.legend()
        plt.show()
    else:
        with open("10.txt", "r") as file:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt = eval(file.read())
        mode="precise"
        if mode=="rough":
            plot_phase_diagram_rough(xt)
        else:
            signallist=[2,0,0,0]
            plot_phase_diagram(xt,signallist,0.05,1,'test.png','g')
        
