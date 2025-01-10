#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.interpolate import splrep

class plotting:
    def __init__(self):
        print("move everything to class later")
        self.phblist=[]

    def plot_phase_diagram(self,phase_boundary_list,signallist,dT,Tstart,filename="phasediagram.png",pointfile="point.txt",color='r'):
        print("plot phase diagram using data obtained from FYL-CVM")
        listlen=len(phase_boundary_list)
        signal_list=signallist.copy()         #protect signallist from modification
        print(phase_boundary_list[1][2][-1])
        print(phase_boundary_list[3][2][-1])
        eutecticend=[[0.4],[0.45],[589]]
        slopelist=[[0,0]]*listlen
        for i in range(listlen):                #search invariant point first
            if signal_list[i]==0:               #a phb that ends
                for j in range(i+1,listlen):
                    if np.abs(phase_boundary_list[i][2][-1]-phase_boundary_list[j][2][-1])<1.0e-6:     #is invariant point
                        slopelist[i],slopelist[j]=self.solve_invariant(phase_boundary_list[i],phase_boundary_list[j],dT,color)
                        signal_list[i]=signal_list[j]=4
        for i in range(listlen):                #finish all endpoint
            if signal_list[i]==0:   
                slopelist[i]=self.solve_endpoint(phase_boundary_list[i],0.5,dT,color)
        for i in range(listlen):
            for j in range(i+1,listlen):       #search for eutectic start
                if phase_boundary_list[i][2][0]==phase_boundary_list[j][2][0] and phase_boundary_list[i][2][0]!=Tstart:
                    if phase_boundary_list[i][0][0]<phase_boundary_list[j][0][0]:
                        eutecticend=self.solve_eutectic(phase_boundary_list[i],phase_boundary_list[j],slopelist[i],slopelist[j],dT,color)
                    else:
                        #print("here is the case")
                        eutecticend=self.solve_eutectic(phase_boundary_list[j],phase_boundary_list[i],slopelist[j],slopelist[i],dT,color)
        for i in range(listlen):
            if signal_list[i]==2:  
                self.finish_eutectic(phase_boundary_list[i],eutecticend,color)
        plt.xlim(0.1, 0.5)
        plt.ylim(1.0, 2.4)
        plt.xlabel("composition")
        plt.ylabel("normalized T")
        if filename:
            plt.savefig(filename)
        plt.show()
        f=open(pointfile,"w")
        print(self.phblist,file=f)
        f.close()
        utility.replace_word_in_file(pointfile,"array","np.array")
        plt.clf()
    
    def plot_phase_boundary_v1(self,xt,color='r',pointfile="point.txt"):
        sig=check_positive_v1(xt[0])
        if sig==1:
            phb=CubicSpline(xt[0],xt[1],bc_type=((2, 0.0), (1, 0.0)))
            xspace=np.linspace(xt[0][0],xt[0][-1],100,endpoint=True)
            plt.plot(xspace,phb(xspace),color)
            self.phblist.append(np.vstack((xspace,phb(xspace))))
            return phb(xt[0][0],1)
        elif sig==-1:
            phb=CubicSpline(1-xt[0],xt[1],bc_type=((2, 0.0), (1, 0.0)))
            xspace=np.linspace(xt[0][0],xt[0][-1],100,endpoint=True)
            plt.plot(xspace,phb(1-xspace),color)
            self.phblist.append(np.vstack((xspace,phb(1-xspace))))
            return -phb(1-xt[0][0],1)
        elif sig==0:                     #plot in yx space
            w=np.ones(len(xt[0]))
            w[-1]=10
            phb=UnivariateSpline(xt[1],xt[0],w=w)
            Tspace=np.linspace(xt[1][0],xt[1][-1],100,endpoint=True)
            plt.plot(phb(Tspace),Tspace,color)
            self.phblist.append(np.vstack((phb(Tspace),Tspace)))
            return 0

    def solve_invariant(self,x_t1,x_t2,dT,color='r'):
        k=((x_t1[0][-1]-x_t2[1][-1])/(x_t1[0][-2]-x_t2[1][-2]))
        T=x_t1[2][-1]+dT*k**(1.35-0.5*k)
        #T=1.8525
        r=0.5*(np.abs((x_t1[0][-1]-x_t1[0][-2])/(x_t2[1][-1]-x_t2[1][-2]))+np.abs((x_t1[1][-1]-x_t1[1][-2])/(x_t2[0][-1]-x_t2[0][-2])))
        x=x_t1[1][-1]+(r/(r+1))*(x_t2[0][-1]-x_t1[1][-1])
        print("invariant point estimation is x={:.3f} T={:.4f}".format(x,T))
        #with open("inv.txt", "a") as file:
        #    file.write(f"{x} {T}\n")
        point=np.array([[x],[x],[T]])                #solve x based on T later, might need to create one FYLCVM class
        slope1=self.plot_invariant(x_t1,point,color)
        slope2=self.plot_invariant(x_t2,point,color)
        return slope1,slope2

    def plot_invariant(self,xt,point,color='r'):
        xt=np.append(xt,point,axis=1)
        slope1=self.plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
        slope2=self.plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)
        return np.array([slope1,slope2])
        #plt.show()

    def solve_endpoint(self,xt,point,dT=10,color='r'):
        Tend=xt[2][-1]+dT*((xt[0][-1]-point)/(xt[0][-2]-point))
        print("Tend is "+str(Tend))
        xt=np.append(xt,np.array([[point],[point],[Tend]]),axis=1)
        slope1=self.plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
        slope2=self.plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)
        return np.array([slope1,slope2])

    def solve_eutectic(self,x_t1,x_t2,slope1,slope2,dT,color='r'):
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
            self.phblist.append(np.vstack((xspace10,Tspace)))
            xspace11=[x_t1[1][0],x1other]
            plt.plot(xspace11,Tspace,color)
            self.phblist.append(np.vstack((xspace11,Tspace)))
            xspace20=[x_t2[0][0],x2other]
            plt.plot(xspace20,Tspace,color)
            self.phblist.append(np.vstack((xspace20,Tspace)))
            xspace21=[x_t2[1][0],x]
            plt.plot(xspace21,Tspace,color)
            self.phblist.append(np.vstack((xspace21,Tspace)))
            xline=[x1other,x,x2other]
            Tline=[T,T,T]
            plt.plot(xline,Tline,color)
            self.phblist.append(np.vstack((xline,Tline)))
            print("T is "+str(T))
            print(np.array([[x2other],[x1other],[T]]))
            return np.array([[x2other],[x1other],[T]])

    def finish_eutectic(self,xt,point,color='r'):
        if np.abs(xt[2][-1]-point[2][0])>1.0e-6:
            print("unequal")
            xt=np.append(xt,point,axis=1)
        else:
            #drop the last element and replace
            xt=xt[:, :-1]
            xt=np.append(xt,point,axis=1)
        slope1=self.plot_phase_boundary_v1(np.vstack((xt[0],xt[2])),color)
        slope2=self.plot_phase_boundary_v1(np.vstack((xt[1],xt[2])),color)

def plot_phase_diagram(phase_boundary_list,signallist,dT,Tstart,filename="phasediagram.png",color='r'):
    print("plot phase diagram using data obtained from FYL-CVM")
    listlen=len(phase_boundary_list)
    signal_list=signallist.copy()         #protect signallist from modification
    print(phase_boundary_list[1][2][-1])
    print(phase_boundary_list[3][2][-1])
    eutecticend=[[0.4],[0.45],[589]]
    slopelist=[[0,0]]*listlen
    for i in range(listlen):                #search invariant point first
        if signal_list[i]==0:               #a phb that ends
            for j in range(i+1,listlen):
                if np.abs(phase_boundary_list[i][2][-1]-phase_boundary_list[j][2][-1])<1.0e-6:     #is invariant point
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
    plt.xlim(0.1, 0.5)
    plt.ylim(1, 2.43)
    plt.xlabel("composition")
    plt.ylabel("normalized T")
    if filename:
        plt.savefig(filename)
    #plt.show()
    #utility.replace_word_in_file("point.txt","array","np.array")
    #plt.clf()

def solve_invariant(x_t1,x_t2,dT,color='r'):
   k=((x_t1[0][-1]-x_t2[1][-1])/(x_t1[0][-2]-x_t2[1][-2]))
   T=x_t1[2][-1]+dT*k**(1.4-0.5*k)
   #T=1.8525
   r=0.5*(np.abs((x_t1[0][-1]-x_t1[0][-2])/(x_t2[1][-1]-x_t2[1][-2]))+np.abs((x_t1[1][-1]-x_t1[1][-2])/(x_t2[0][-1]-x_t2[0][-2])))
   x=x_t1[1][-1]+(r/(r+1))*(x_t2[0][-1]-x_t1[1][-1])
   print("invariant point estimation is x={:.3f} T={:.4f}".format(x,T))
   #with open("inv.txt", "a") as file:
   #    file.write(f"{x} {T}\n")
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
        with open("11.txt", "r") as file1:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt = eval(file1.read())
        plot_phase_diagram(xt,signallist,0.05,1,False,'g')
        print(signallist)
        with open("10.txt", "r") as file2:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt2 = eval(file2.read())
        plot_phase_diagram(xt2,signallist,0.05,1,False,'r')
        with open("09.txt", "r") as file3:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt3 = eval(file3.read())
        plot_phase_diagram(xt3,signallist,0.05,1,False,'b')
        plt.xlim(0.1, 0.5)
        plt.ylim(1, 2.4)
        plt.xlabel("x")
        plt.ylabel("normalized T")
        plt.plot(0,0,color='r',label="r=1")
        plt.plot(0,0,color='g',label="r=1.1")
        plt.plot(0,0,color='b',label="r=0.9")
        plt.legend()
        plt.show()
        '''
        signallist=np.array([2,0,0,0])
        colors=['r','orange','y','g','b']
        print("shell mode, input whatever you want")
        for i in range(5):
            p=1+0.02*i
            filename="local"+str(p)+".txt"
            #figurename="local"+str(i)+".png"
            figurename="zarbage.png"
            with open(filename, "r") as file:
                xt = eval(file.read())
                plot_phase_diagram(xt,signallist,0.05,1,figurename,colors[i])

        red_patch = mpatches.Patch(color='red', label='p=1')
        blue_patch = mpatches.Patch(color='blue', label='p=1.08')
        plt.legend(handles=[red_patch, blue_patch])
        plt.show()'''
        
    else:
        signallist=[2,0,0,0]
        with open("p1.05.txt", "r") as file:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt = eval(file.read())
        mode="class"
        if mode=="rough":
            plot_phase_diagram(xt)
        elif mode=="class":
            myclass=plotting()
            myclass.plot_phase_diagram(xt,signallist,0.05,1,'placeholder.png','point1.05.txt','g')
        else:
            plot_phase_diagram(xt,signallist,0.05,1,'placeholder.png','g')
        
