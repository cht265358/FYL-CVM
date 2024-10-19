#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
from scipy.optimize import minimize

#This is the utility function file include all the "dirty" funtions

def get_cluster_prob(high_d_matrix,position_type_matrix):     #get cluster probablity from high dimensional matrix
    if np.abs(np.sum(high_d_matrix)-1.0)>0.0001:
        print("error, total probability is not equal to 1 but %f",np.sum(high_d_matrix))
        return 0
    shape=np.shape(high_d_matrix)
    prob=0
    if np.shape(position_type_matrix)[1]==1:     #compute point cluster probability
        location=position_type_matrix[0][0]
        atom_type=position_type_matrix[1][0]
        for idx in itertools.product(*[range(s) for s in shape]):
            if idx[int(location)]==atom_type:
                prob+=high_d_matrix[idx]
    if np.shape(position_type_matrix)[1]==2:     #compute binary cluster probability
        location=position_type_matrix[0]
        atom_type=position_type_matrix[1]
        for idx in itertools.product(*[range(s) for s in shape]):
            if idx[int(location[0])]==atom_type[0] and idx[int(location[1])]==atom_type[1]:
                prob+=high_d_matrix[idx]
    return prob

class userinput:             #read in user input

    def __init__(self,usrinput):
        inputdictionary={
            "-t":[False],
            "main.py":[True],
            "-c":[2,"int"],
            "-maxsize":[4,"int"],
            "-cs":["FCC","string"],
            "-basic":["basic.in","string"],
            "-point":[False],
            "-vib":[1,"float"],
            "-e":[0,"float"],
            "-h":[False],
            "-disp":[False],
            "-name":["phasediagram.png","string"],
            "-dmu":[0.05,"float"],
            "-unit":["hypo","string"],
            "-control":["control.txt","string"],
        }
        self.inputdictionary=inputdictionary
        self.usrinput=usrinput
    
    def read_input(self):
        for keyword in self.usrinput:
            #print(keyword)
            if detectequal(keyword):
                inputtuple=keyword.split("=")
                if self.inputdictionary[inputtuple[0]][1]=="int":
                    self.inputdictionary.update({inputtuple[0]:[int(inputtuple[1]),"int"]})
                elif self.inputdictionary[inputtuple[0]][1]=="float":
                    self.inputdictionary.update({inputtuple[0]:[float(inputtuple[1]),"float"]})
                elif self.inputdictionary[inputtuple[0]][1]=="string":
                    self.inputdictionary.update({inputtuple[0]:[inputtuple[1],"string"]})
            else:
                self.inputdictionary.update({keyword:True})

def detectequal(astring):
    for element in astring:
        if element=="=":
            return 1
    return 0

def read_control_file(controlfilename):
    control_dictionary={
        "dT":0.05,             
        "dmurough":0.1,
        "dmuprecise":0.02
    }
    return control_dictionary

def read_control_variable(controlfilename,unit):
    if os.path.isfile(controlfilename):
        control_dictionary=read_control_file(controlfilename)
        return control_dictionary
    else:          #use default value based on unit
        if unit=="hypo":
            control_dictionary={
                "dT":0.05,             
                "dmurough":0.1,
                "dmuprecise":0.02,
                "dmuscan":0.2,
                "Tmax":2.5,
                "Tstart":1.0,
                "R":1.0
            }
            return control_dictionary
        elif unit=="kj":
            control_dictionary={              #adjust later for real unit system
                "dT":10.0,             
                "dmurough":0.1,
                "dmuprecise":0.02,
                "dmuscan":0.2,
                "Tmax":1000.0,
                "Tstart":500.0,           #change to 300 later
                "R":8.3144621e-3
            }
            return control_dictionary

def displayhelptext():
    print("Hello world, this is help information")
    sys.exit()
    
if __name__ == '__main__':
    arr = np.ones([2,2,2,2])
    arr = arr/np.sum(arr)
    print(np.sum(arr))
    position_type_matrix=np.array([[2,3],[1,1]],ndmin=2)
    print(get_cluster_prob(arr,position_type_matrix))


'''

    def trace_phase_boundary(self,Tmin,Tmax,steplength,mustart=0):        #trace the phase boundary between two phases, will terminate when new phase shows up
        x1mat=np.zeros(round((Tmax-Tmin)/steplength)+1)
        print(x1mat)
        x2mat=np.zeros(round((Tmax-Tmin)/steplength)+1)
        Tspace=np.linspace(Tmin,Tmax,round((Tmax-Tmin)/steplength+1),endpoint=True)
        T=Tmin
        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,0.1,500)      #set to 7 to speed up calculation
        print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2))
        #use finer grid to locate two phase eqm
        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart-0.01,0.01,12)
        print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2))
        print("sitevar1 is "+str(result1.x)+" sitevar2 is "+str(result2.x))  
        lowphase=self.identify_phase(result1.x)             #phase with low potential value, these two value won't change when T increases
        highphase=self.identify_phase(result2.x)             #phase with high potential value
        x1mat[0]=x1
        x2mat[0]=x2
        count=1
        while (T<=Tmax):
            T+=steplength
            mustart+=self.compute_dmu()
            result=self.optimize_grand_potential(T,mustart,"BFGS")
            currentphase=self.identify_phase(result.x)
            print("current phase is "+currentphase)
            if currentphase==lowphase:
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,0.01,100)
            elif currentphase==highphase:
                #search in reverse so that result is also in reverse
                (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-0.01,100) 
            #search with higher precision, we can add a loop later
            #(mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart-0.001,0.001,100)
            print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2)+" at T="+str(T))
            print("sitevar1 is "+str(result1.x)+" sitevar2 is "+str(result2.x))     
            x1mat[count]=x1
            x2mat[count]=x2  
            count+=1
        plt.plot(x1mat,Tspace)
        plt.plot(x2mat,Tspace)
        #plt.show()   
    '''