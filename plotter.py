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

def plot_phase_diagram(phase_boundary_list,mu_list,signal_list,dT,Tstart,data_type="phaseboundary"):
    print("plot phase diagram using data obtained from FYL-CVM")
    listlen=len(phase_boundary_list)
    eutecticend=[[0.4],[0.45],[589]]
    slopelist=[[0,0]]*listlen
    for i in range(listlen):                #search invariant point first
        if signal_list[i]==0:               #a phb that ends
            for j in range(i+1,listlen):
                if phase_boundary_list[i][2][-1]==phase_boundary_list[j][2][-1]:     #is invariant point
                    slopelist[i],slopelist[j]=solve_invariant(phase_boundary_list[i],phase_boundary_list[j],mu_list[i],mu_list[j])
                    signal_list[i]=signal_list[j]=4
    for i in range(listlen):                #finish all endpoint
        if signal_list[i]==0:   
            slopelist[i]=solve_endpoint(phase_boundary_list[i],0.5,dT)
    for i in range(listlen):
        for j in range(i+1,listlen):       #search for eutectic start
            if phase_boundary_list[i][2][0]==phase_boundary_list[j][2][0] and phase_boundary_list[i][2][0]!=Tstart:
                if phase_boundary_list[i][0][0]<phase_boundary_list[j][0][0]:
                    eutecticend=solve_eutectic(phase_boundary_list[i],phase_boundary_list[j],slopelist[i],slopelist[j])
                else:
                    #print("here is the case")
                    eutecticend=solve_eutectic(phase_boundary_list[j],phase_boundary_list[i],slopelist[j],slopelist[i])
    for i in range(listlen):
        if signal_list[i]==2:  
            finish_eutectic(phase_boundary_list[i],eutecticend)
    plt.show()

def solve_invariant(x_t1,x_t2,mu_t1,mu_t2):
    #x_t1,flipx1=check_positive(x_t1)
    #x_t2,flipx2=check_positive(x_t2)
    mu,T=find_intersect(mu_t1,mu_t2)
    point=np.array([[0.265],[0.265],[2.03]])                #solve x based on T later, might need to create one FYLCVM class
    slope1=plot_invariant(x_t1,point)
    slope2=plot_invariant(x_t2,point)
    return slope1,slope2
    #plt.show()

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
    print("solve end point now")
    xt,flip=check_positive(xt)
    if flip==1:
        #phb2=UnivariateSpline(xt[1],xt[2])              #change to other spline if necessary
        #Tend=phb2(point)
        Tend=xt[2][-1]+0.5*dT
        xt=np.append(xt,np.array([[point],[point],[Tend]]),axis=1)
        phb1=CubicSpline(xt[0],xt[2],bc_type=((2, 0.0), (1, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((2, 0.0), (1, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')
        return np.array([phb1(xt[0][0],1),phb2(xt[1][0],1)])
    elif flip==-1:
        Tend=xt[2][0]+0.5*dT
        xt=np.append(np.array([[point],[point],[Tend]]),xt,axis=1)
        phb1=CubicSpline(xt[0],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        phb2=CubicSpline(xt[1],xt[2],bc_type=((1, 0.0), (2, 0.0)))
        x1=np.linspace(xt[0][0],xt[0][-1],1000,endpoint=True)
        x2=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
        plt.plot(x1,phb1(x1),'r')
        plt.plot(x2,phb2(x2),'r')
        return np.array([phb1(xt[0][-1],1),phb2(xt[1][-1],1)])

def solve_eutectic(x_t1,x_t2,slope1,slope2):
    if x_t1[0][0]<x_t2[0][0]:                   #make sure xt1 on the left
        Tend=x_t1[2][0]
        x=(slope1[0]*x_t1[0][0]-slope2[1]*x_t2[1][0])/(slope1[0]-slope2[1])
        T=x_t1[2][0]+(x-x_t1[0][0])*slope1[0]
        x1other=x_t1[1][0]+(T-x_t1[2][0])/slope1[1]
        x2other=x_t2[0][0]+(T-x_t2[2][0])/slope2[0]
        Tspace=[Tend,T]
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
        return np.array([[x2other],[x1other],[T]])

def finish_eutectic(xt,point):
    xt=np.append(xt,point,axis=1)
    #print(xt)
    xright=np.vstack((xt[0],xt[2]))
    print(xright)
    xleft=np.vstack((xt[1],xt[2]))
    xright,flipright=check_positive(xright)
    xleft,flipleft=check_positive(xleft)
    print(xright)
    if flipright==1:
        phbright=CubicSpline(xright[0],xright[1],bc_type=((2, 0.0), (1, 0.0)))
    else:
        phbright=CubicSpline(xright[0],xright[1],bc_type=((1, 0.0), (2, 0.0)))
    if flipleft==1:
        phbleft=CubicSpline(xleft[0],xleft[1],bc_type=((2, 0.0), (1, 0.0)))
    else:
        phbleft=CubicSpline(xleft[0],xleft[1],bc_type=((1, 0.0), (2, 0.0)))
    right=np.linspace(xright[0][0],xright[0][-1],1000,endpoint=True)
    left=np.linspace(xleft[0][0],xleft[0][-1],1000,endpoint=True)
    plt.plot(right,phbright(right),'r')
    plt.plot(left,phbleft(left),'r')

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

    phb_list=[np.array([[0.40910661, 0.40863414, 0.40868302, 0.40854587, 0.4088605 ,
        0.40833552, 0.40883314, 0.40907053, 0.41022444, 0.41102828,
        0.41323821, 0.41620877, 0.42055948],
       [0.3698831 , 0.37120107, 0.37304854, 0.37492779, 0.37726884,
        0.3791905 , 0.3819503 , 0.38468543, 0.38816119, 0.39157075,
        0.39599599, 0.4009904 , 0.40684768],
       [1.        , 1.05      , 1.1       , 1.15      , 1.2       ,
        1.25      , 1.3       , 1.35      , 1.4       , 1.45      ,
        1.5       , 1.55      , 1.6       ]]), np.array([[0.1983263 , 0.19946549, 0.2002097 , 0.20163571, 0.20278955,
        0.20450277, 0.20666342, 0.20847593, 0.21062881, 0.21243494,
        0.21450237, 0.21674449, 0.21910341, 0.22154489, 0.22433734,
        0.22712248, 0.2303955 , 0.23381772, 0.23782346, 0.24316181,
        0.25100075],
       [0.13962165, 0.14096667, 0.14321775, 0.14550723, 0.14783303,
        0.15057211, 0.1537    , 0.15682954, 0.16032082, 0.16380354,
        0.16762857, 0.1717833 , 0.176256  , 0.18103981, 0.18645678,
        0.19216302, 0.1988029 , 0.20603382, 0.21448235, 0.22539103,
        0.24058251],
       [1.        , 1.05      , 1.1       , 1.15      , 1.2       ,
        1.25      , 1.3       , 1.35      , 1.4       , 1.45      ,
        1.5       , 1.55      , 1.6       , 1.65      , 1.7       ,
        1.75      , 1.8       , 1.85      , 1.9       , 1.95      ,
        2.        ]]), np.array([[0.39948941, 0.37929707, 0.36451222, 0.3516783 , 0.33922487,
        0.3265239 , 0.31202252, 0.29356041],
       [0.39836591, 0.37685969, 0.36060341, 0.3463028 , 0.33246253,
        0.31872417, 0.30384841, 0.28685345],
       [1.65      , 1.7       , 1.75      , 1.8       , 1.85      ,
        1.9       , 1.95      , 2.        ]]), np.array([[0.42965505, 0.44027972, 0.45048576, 0.45988982, 0.470064  ,
        0.48139017],
       [0.4151361 , 0.42364311, 0.43273081, 0.44239539, 0.45418061,
        0.46960722],
       [1.65      , 1.7       , 1.75      , 1.8       , 1.85      ,
        1.9       ]])]
    muTlist=[np.array([[7.45, 7.41, 7.35, 7.29, 7.21, 7.15, 7.05, 6.95, 6.81, 6.67, 6.47,
        6.23, 5.93],
       [1.  , 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.5 ,
        1.55, 1.6 ]]), np.array([[22.53, 22.45, 22.33, 22.21, 22.09, 21.95, 21.79, 21.63, 21.45,
        21.27, 21.07, 20.85, 20.61, 20.35, 20.05, 19.73, 19.35, 18.93,
        18.43, 17.77, 16.83],
       [ 1.  ,  1.05,  1.1 ,  1.15,  1.2 ,  1.25,  1.3 ,  1.35,  1.4 ,
         1.45,  1.5 ,  1.55,  1.6 ,  1.65,  1.7 ,  1.75,  1.8 ,  1.85,
         1.9 ,  1.95,  2.  ]]), np.array([[ 6.49,  7.79,  8.75,  9.59, 10.41, 11.25, 12.21, 13.43],
       [ 1.65,  1.7 ,  1.75,  1.8 ,  1.85,  1.9 ,  1.95,  2.  ]]), np.array([[5.47, 4.93, 4.35, 3.73, 2.97, 1.97],
       [1.65, 1.7 , 1.75, 1.8 , 1.85, 1.9 ]])]
    signallist=[2,0,0,0]
    #solve_invariant(testinputx[0],testinputx[1],testinput[0],testinput[1])
    plot_phase_diagram(phb_list,muTlist,signallist,0.05,1)