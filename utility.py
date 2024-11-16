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
            "-shell":[False],
            "main.py":[True],
            "-c":[2,"int"],
            "-maxsize":[4,"int"],
            "-cs":["FCC","string"],
            "-basic":["basic.in","string"],
            "-point":[False],
            "-vib":[1.0,"float"],
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
        self.read_input()
    
    def assign_input(self):
        self.testmode=self.inputdictionary['-t']
        self.shellmode=self.inputdictionary['-shell']
        self.number_of_component=self.inputdictionary['-c'][0]
        self.max_clustersize=self.inputdictionary['-maxsize'][0]
        self.crystal_system=self.inputdictionary['-cs'][0]
        self.basic_cluster_energy,self.local_energy_parameter=read_cluster_energy(self.inputdictionary['-basic'][0])
        self.vibration_parameter=self.inputdictionary['-vib'][0]
        self.displayhelptext=self.inputdictionary['-h']
        self.elastic_parameter=self.inputdictionary['-e'][0]
        self.phasediagram_name=self.inputdictionary['-name'][0]
        #potential_precision=inputs.inputdictionary['-dmu'][0]
        self.unituse=self.inputdictionary['-unit'][0]
        self.controlfilename=self.inputdictionary['-control'][0]
    
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
        self.assign_input()

class phaseboundary:
    def __init__(self,Tstart,mustart,muend,lowphase,highphase,c1,c2,direction=1):
        self.status="open"
        self.lowphase=lowphase
        self.highphase=highphase
        self.mustart=mustart
        self.muend=muend
        self.mustartpoint=0.5*(mustart+muend)
        self.muendpoint=0.0
        self.Tstart=Tstart
        self.Tend=0
        self.x1mat=np.array([c1])
        self.x2mat=np.array([c2])
        self.Tspace=np.array([Tstart]) 
        self.muspace=np.array([0.5*(mustart+muend)])
        self.direction=direction
        print("creating new node between "+lowphase+" and "+highphase+" at T="+str(Tstart))

def compare_phb(phb1,phb2,mutol=0.1,Ttol=0.005):   #phb1 is the phb to be determined and phb2 is the existing phb
    #return 1 means there is overlap with other nodes
    if phb2.status=="exit" or phb2.status=="boundary":
        phb2.muendpoint=phb2.muspace[-1]
    else:
        return 0
    if np.abs(phb1.mustartpoint-phb2.mustartpoint)<mutol and np.abs(phb1.Tstart-phb2.Tstart)<Ttol:
        return 1      
    elif np.abs(phb1.mustartpoint-phb2.muendpoint)<mutol and np.abs(phb1.Tstart-phb2.Tend)<Ttol:
        return 1
    else:
        return 0

def order_phb(phb):
    if phb[0][1]<phb[0][0]:
        return np.flip(phb,axis=0)
    else:
        return phb
    
def read_cluster_energy(filename):
    txt=np.loadtxt(filename)
    if txt.ndim==1:
        return txt,np.ones(np.shape(txt)[0])
    else:
        return txt[0],txt[1]

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
                "dmurough":0.2,
                "dmuprecise":0.01,
                "dmuscan":1.0,
                "Tmax":3.5,
                "Tstart":1.0,
                "R":1.0,
                "mutol":0.0101,
                "mumin":0.0,
                "mumax":25.0
            }
            return control_dictionary
        elif unit=="kj":
            control_dictionary={              #adjust later for real unit system
                "dT":25.0,             
                "dmurough":0.5,
                "dmuprecise":0.05,
                "dmuscan":0.2,
                "Tmax":900.0,
                "Tstart":550.0,           #change to 300 later
                "R":8.3144621e-3,
                "mutol":0.0201,
                "mumin":-75.0,
                "mumax":75.0
            }
            return control_dictionary

def displayhelptext():
    print("Hello world, this is help information")
    sys.exit()

def map_variable_to_phase(initialvalue,phase,n):
    sitepotential_input=np.zeros(n)
    if phase=="L12":
        sitepotential_input[0]=sitepotential_input[1]=sitepotential_input[2]=initialvalue[0]
        sitepotential_input[3]=initialvalue[1]
    elif phase=="L10":
        sitepotential_input[0]=sitepotential_input[1]=initialvalue[0]
        sitepotential_input[2]=sitepotential_input[3]=initialvalue[1]        
    return sitepotential_input


def replace_word_in_file(file_path, old_word, new_word):
   try:
       # Open the file in read mode to get its content
       with open(file_path, 'r') as file:
           content = file.read()

       # Replace all instances of the old word with the new word
       updated_content = content.replace(old_word, new_word)

       # Open the file in write mode to update it with the new content
       with open(file_path, 'w') as file:
           file.write(updated_content)
               
   except FileNotFoundError:
       print(f"Error: The file '{file_path}' does not exist.")
   except IOError:
       print("An error occurred while reading or writing the file.")

def compute_dmu(muspace):
    if len(muspace)<=1:
        return 0
    else:
        return muspace[-1]-muspace[-2]

if __name__ == '__main__':
    a,b=read_cluster_energy("basic.in")
    print(b)


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