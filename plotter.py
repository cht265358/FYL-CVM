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

def plot_phase_diagram(phase_boundary_list,signal_list,dT,data_type="phaseboundary"):
    print("plot phase diagram using data obtained from FYL-CVM")
    listlen=len(phase_boundary_list)
    
    for i in range(listlen):
        for j in range(i+1,listlen):
            if signal_list[i]!=4 and signal_list[j]!=4:
                if phase_boundary_list[i][2][-1]==phase_boundary_list[j][2][-1]:    #same end point
                    #if phase_boundary_list[i][0][0]>
                    phase_boundary_list[j][[0,1]]=phase_boundary_list[j][[1,0]]
                    list_use=np.append(phase_boundary_list[i],phase_boundary_list[j],axis=1)
                    if list_use[0][0]>list_use[0][1]:
                        list_use=np.fliplr(list_use) 
                    phb_right=UnivariateSpline(list_use[0],list_use[2],s=0.5)
                    phb_left=UnivariateSpline(list_use[1],list_use[2],s=0.5)
                    x1space=np.linspace(list_use[0][0],list_use[0][-1],1000)
                    x2space=np.linspace(list_use[1][0],list_use[1][-1],1000)
                    plt.plot(x1space,phb_right(x1space),'r')
                    plt.plot(x2space,phb_left(x2space),'r')
                    signal_list[i]=signal_list[j]=4
                    #if phase_boundary_list[i][1][-1]>phase_boundary_list[j][0][-1]:
                    #avg=0.5*(phase_boundary_list[i][1][-1]+phase_boundary_list[j][0][-1])
                    #phase_boundary_list[i]=np.append(phase_boundary_list[i],np.array([[avg],[avg],[phase_boundary_list[i][2][-1]+0.5*dT]]),axis=1)
                    #phase_boundary_list[i]=np.append(phase_boundary_list[i],phase_boundary_list[j],axis=1)
                    #phase_boundary_list[j]=np.append(phase_boundary_list[j],np.array([[avg],[avg],[phase_boundary_list[i][2][-1]]]),axis=1)
    '''
    print(phase_boundary_list)
    for i in range(listlen):
        list_use=phase_boundary_list[i]
        if len(list_use[0])>3:
            phb_right=UnivariateSpline(list_use[0],list_use[2],s=0.5)
            phb_left=UnivariateSpline(list_use[1],list_use[2],s=0.5)
            x1space=np.linspace(list_use[0][0],list_use[0][-1],1000)
            x2space=np.linspace(list_use[0][0],list_use[0][-1],1000)
            plt.plot(x1space,phb_right(x1space),'r')
            plt.plot(x2space,phb_left(x2space),'r')
        else:
            plt.plot(list_use[0],list_use[2],'r')
            plt.plot(list_use[1],list_use[2],'r')
    plt.show()
    '''

def solve_invariant(x_t1,x_t2,mu_t1,mu_t2):
    #x_t1,flipx1=check_positive(x_t1)
    #x_t2,flipx2=check_positive(x_t2)
    mu,T=find_intersect(mu_t1,mu_t2)
    point=np.array([[0.732],[0.732],[763.0606331864144]])
    plot_invariant(x_t1,point)
    plot_invariant(x_t2,point)
    plt.show()

def plot_invariant(xt,point):
    xt,flip=check_positive(xt)
    if flip==1:
        xt=np.append(xt,point,axis=1)
        phb1=CubicSpline(xt[0],xt[2],bc_type=((2, 0.0), (1, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((2, 0.0), (1, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')
        return np.array([phb1(xt[0][0],1),phb2(xt[1][0],1)])
        #plt.show()
    elif flip==-1:
        xt=np.append(point,xt,axis=1)
        phb1=CubicSpline(xt[0],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')
        return np.array([phb1(xt[0][-1],1),phb2(xt[1][-1],1)])
        #plt.show()

def solve_endpoint(xt,point,dT=10):
    xt,flip=check_positive(xt)
    if flip==1:
        #phb2=UnivariateSpline(xt[1],xt[2])              #change to other spline if necessary
        #Tend=phb2(point)
        Tend=xt[2][-1]+0.5*dT
        xt=np.append(xt,np.array([[point],[point],[Tend]]))
        phb1=CubicSpline(xt[0],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')
    elif flip==-1:
        Tend=xt[2][0]+0.5*dT
        xt=np.append(np.array([[point],[point],[Tend]]),xt)
        phb1=CubicSpline(xt[0],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')

def solve_eutectic(x_t1,x_t2,slope1,slope2):
    if x_t1[0][0]<x_t2[0][0]:                   #make sure xt1 on the left
        Tend=x_t1[2][0]
        x=(slope1[0]*x_t1[0][0]-slope2[1]*x_t2[1][0])/(slope1[0]-slope2[1])
        T=x_t1[2][0]+(x-x_t1[0][0])*slope1[0]
        x1other=x_t1[1][0]+(T-x_t1[2][0])/slope1[1]
        x2other=x_t2[0][0]+(T-x_t2[2][0])/slope2[0]
        Tspace=[T,Tend]
        xspace10=[x_t1[0][0],x]
        plt.plot(xspace10,Tspace,'g')
        xspace11=[x_t1[1][0],x1other]
        plt.plot(xspace11,Tspace,'g')
        xspace20=[x_t2[0][0],x2other]
        plt.plot(xspace20,Tspace,'g')
        xspace21=[x_t2[1][0],x]
        plt.plot(xspace21,Tspace,'g')
        xline=[x1other,x,x2other]
        Tline=[T,T,T]
        plt.plot(xline,Tline,'g')
        return np.array([[x2other],[x1other],T])

def finish_eutectic(xt,point):
    xt=np.append()

def find_intersect(mu_t1,mu_t2):
    mu_t1,flipmu1=check_positive(mu_t1)
    mu_t2,flipmu2=check_positive(mu_t2)
    fit1=UnivariateSpline(mu_t1[0],mu_t1[1])
    fit2=UnivariateSpline(mu_t2[0],mu_t2[1])
    Fmin=10
    mubest=0
    Tbest=300
    if flipmu1+flipmu2==2:    #two positive line
        ls=np.linspace(np.max(mu_t1[0][-1],mu_t2[0][-1]),np.max(mu_t1[0][-1],mu_t2[0][-1])+10,1001,endpoint=True)
    elif flipmu1==1 and flipmu2==-1:
        ls=np.linspace(mu_t1[0][-1],mu_t2[0][0],500)
    elif flipmu1==-1 and flipmu2==1:
        ls=np.linspace(mu_t2[0][-1],mu_t1[0][0],500)
    else:
        ls=np.linspace(np.min(mu_t1[0][0],mu_t2[0][0]),np.min(mu_t1[0][0],mu_t2[0][0])-10,1001,endpoint=True)
    for i in ls:
        F=np.abs(fit1(i)-fit2(i))
        if F<Fmin:
            Fmin=F
            mubest=i
            Tbest=0.5*(fit1(i)+fit2(i))
    return mubest,Tbest
    


def check_positive(x):
    flip=1
    if x[0][0]>x[0][-1]:
        x=np.fliplr(x)
        flip=-1
    return x,flip

if __name__ == '__main__':

    testinput=[np.array([[-53.925, -53.075, -52.125, -51.075, -49.775, -48.125, -45.475],
       [700.   , 710.   , 720.   , 730.   , 740.   , 750.   , 760.   ]]), np.array([[-28.625, -29.975, -31.375, -32.925, -34.625, -36.725, -39.775],
       [700.   , 710.   , 720.   , 730.   , 740.   , 750.   , 760.   ]])]
    testinputx=[np.array([[7.92327224e-01, 7.87730420e-01, 7.82644248e-01, 7.77074867e-01,
        7.70269436e-01, 7.61737463e-01, 7.48255546e-01],
       [7.65459901e-01, 7.63260499e-01, 7.60816139e-01, 7.58133191e-01,
        7.54747779e-01, 7.50343537e-01, 7.42817140e-01],
       [7.00000000e+02, 7.10000000e+02, 7.20000000e+02, 7.30000000e+02,
        7.40000000e+02, 7.50000000e+02, 7.60000000e+02]]), np.array([[6.73028602e-01, 6.80105253e-01, 6.87207277e-01, 6.94826942e-01,
        7.02750338e-01, 7.11963980e-01, 7.24047493e-01],
       [6.65327222e-01, 6.71902798e-01, 6.78710927e-01, 6.86245660e-01,
        6.94509752e-01, 7.04740648e-01, 7.19668575e-01],
       [7.00000000e+02, 7.10000000e+02, 7.20000000e+02, 7.30000000e+02,
        7.40000000e+02, 7.50000000e+02, 7.60000000e+02]])]
    signallist=[0,0]
    solve_invariant(testinputx[0],testinputx[1],testinput[0],testinput[1])
    #plot_phase_diagram(testinput,signallist,5.0)