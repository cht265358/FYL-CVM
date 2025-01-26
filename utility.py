#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
from scipy.optimize import minimize

#This is the utility function file include all the "dirty" funtions

def delta(i,j):
    return 1 if i == j else 0

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
    elif np.shape(position_type_matrix)[1]==2:     #compute binary cluster probability
        location=position_type_matrix[0]
        atom_type=position_type_matrix[1]
        for idx in itertools.product(*[range(s) for s in shape]):
            if idx[int(location[0])]==atom_type[0] and idx[int(location[1])]==atom_type[1]:
                prob+=high_d_matrix[idx]
    elif np.shape(position_type_matrix)[1]==3:     #compute binary cluster probability
        location=position_type_matrix[0]
        atom_type=position_type_matrix[1]
        for idx in itertools.product(*[range(s) for s in shape]):
            if idx[int(location[0])]==atom_type[0] and idx[int(location[1])]==atom_type[1] and idx[int(location[2])]==atom_type[2]:
                prob+=high_d_matrix[idx]
    return prob

def compute_binary_cluster_energy(prob,lattice="BCC"):
    if lattice=="BCC":
        para=np.array([1.5,1])
    elif lattice=="FCC":
        para=np.array([-1,-1])
    energy=0.0
    ab00 = prob[0][0][0][0] + prob[0][0][0][1] + prob[0][0][1][0] + prob[0][0][1][1]
    ab01 = prob[0][1][0][0] + prob[0][1][0][1] + prob[0][1][1][0] + prob[0][1][1][1]
    ab10 = prob[1][0][0][0] + prob[1][0][0][1] + prob[1][0][1][0] + prob[1][0][1][1]
    ab11 = prob[1][1][0][0] + prob[1][1][0][1] + prob[1][1][1][0] + prob[1][1][1][1]
    ac00 = prob[0][0][0][0] + prob[0][0][0][1] + prob[0][1][0][0] + prob[0][1][0][1]
    ac01 = prob[0][0][1][0] + prob[0][0][1][1] + prob[0][1][1][0] + prob[0][1][1][1]
    ac10 = prob[1][0][0][0] + prob[1][0][0][1] + prob[1][1][0][0] + prob[1][1][0][1]
    ac11 = prob[1][0][1][0] + prob[1][0][1][1] + prob[1][1][1][0] + prob[1][1][1][1]
    ad00 = prob[0][0][0][0] + prob[0][0][1][0] + prob[0][1][0][0] + prob[0][1][1][0]
    ad01 = prob[0][0][0][1] + prob[0][0][1][1] + prob[0][1][0][1] + prob[0][1][1][1]
    ad10 = prob[1][0][0][0] + prob[1][0][1][0] + prob[1][1][0][0] + prob[1][1][1][0]
    ad11 = prob[1][0][0][1] + prob[1][0][1][1] + prob[1][1][0][1] + prob[1][1][1][1]
    bc00 = prob[0][0][0][0] + prob[0][0][0][1] + prob[1][0][0][0] + prob[1][0][0][1]
    bc01 = prob[0][0][1][0] + prob[0][0][1][1] + prob[1][0][1][0] + prob[1][0][1][1]
    bc10 = prob[0][1][0][0] + prob[0][1][0][1] + prob[1][1][0][0] + prob[1][1][0][1]
    bc11 = prob[0][1][1][0] + prob[0][1][1][1] + prob[1][1][1][0] + prob[1][1][1][1]
    bd00 = prob[0][0][0][0] + prob[0][0][1][0] + prob[1][0][0][0] + prob[1][0][1][0]
    bd01 = prob[0][0][0][1] + prob[0][0][1][1] + prob[1][0][0][1] + prob[1][0][1][1]
    bd10 = prob[0][1][0][0] + prob[0][1][1][0] + prob[1][1][0][0] + prob[1][1][1][0]
    bd11 = prob[0][1][0][1] + prob[0][1][1][1] + prob[1][1][0][1] + prob[1][1][1][1]
    cd00 = prob[0][0][0][0] + prob[0][1][0][0] + prob[1][0][0][0] + prob[1][1][0][0]
    cd01 = prob[0][0][0][1] + prob[0][1][0][1] + prob[1][0][0][1] + prob[1][1][0][1]
    cd10 = prob[0][0][1][0] + prob[0][1][1][0] + prob[1][0][1][0] + prob[1][1][1][0]
    cd11 = prob[0][0][1][1] + prob[0][1][1][1] + prob[1][0][1][1] + prob[1][1][1][1]

    result1 = (ab00 * np.log(ab00) +ab01 * np.log(ab01) +ab10 * np.log(ab10) +ab11 * np.log(ab11) +
    cd00 * np.log(cd00) +cd01 * np.log(cd01) +cd10 * np.log(cd10) +cd11 * np.log(cd11))
    result2 = (ac00 * np.log(ac00) +ac01 * np.log(ac01) +ac10 * np.log(ac10) +ac11 * np.log(ac11) +
    ad00 * np.log(ad00) +ad01 * np.log(ad01) +ad10 * np.log(ad10) +ad11 * np.log(ad11) +
    bc00 * np.log(bc00) +bc01 * np.log(bc01) +bc10 * np.log(bc10) +bc11 * np.log(bc11) +
    bd00 * np.log(bd00) +bd01 * np.log(bd01) +bd10 * np.log(bd10) +bd11 * np.log(bd11))
    return para[0]*result1+para[1]*result2

def compute_ternary_cluster_energy(prob,lattice="BCC"):
    if lattice=="BCC":
        para=-3.0
    else:
        para=0.0
    abc000 = prob[0][0][0][0] + prob[0][0][0][1]
    abc001 = prob[0][0][1][0] + prob[0][0][1][1]
    abc010 = prob[0][1][0][0] + prob[0][1][0][1]
    abc011 = prob[0][1][1][0] + prob[0][1][1][1]
    abc100 = prob[1][0][0][0] + prob[1][0][0][1]
    abc101 = prob[1][0][1][0] + prob[1][0][1][1]
    abc110 = prob[1][1][0][0] + prob[1][1][0][1]
    abc111 = prob[1][1][1][0] + prob[1][1][1][1]

    abd000 = prob[0][0][0][0] + prob[0][0][1][0]
    abd001 = prob[0][0][0][1] + prob[0][0][1][1]
    abd010 = prob[0][1][0][0] + prob[0][1][1][0]
    abd011 = prob[0][1][0][1] + prob[0][1][1][1]
    abd100 = prob[1][0][0][0] + prob[1][0][1][0]
    abd101 = prob[1][0][0][1] + prob[1][0][1][1]
    abd110 = prob[1][1][0][0] + prob[1][1][1][0]
    abd111 = prob[1][1][0][1] + prob[1][1][1][1]

    acd000 = prob[0][0][0][0] + prob[0][1][0][0]
    acd001 = prob[0][0][0][1] + prob[0][1][0][1]
    acd010 = prob[0][0][1][0] + prob[0][1][1][0]
    acd011 = prob[0][0][1][1] + prob[0][1][1][1]
    acd100 = prob[1][0][0][0] + prob[1][1][0][0]
    acd101 = prob[1][0][0][1] + prob[1][1][0][1]
    acd110 = prob[1][0][1][0] + prob[1][1][1][0]
    acd111 = prob[1][0][1][1] + prob[1][1][1][1]

    bcd000 = prob[0][0][0][0] + prob[1][0][0][0]
    bcd001 = prob[0][0][0][1] + prob[1][0][0][1]
    bcd010 = prob[0][0][1][0] + prob[1][0][1][0]
    bcd011 = prob[0][0][1][1] + prob[1][0][1][1]
    bcd100 = prob[0][1][0][0] + prob[1][1][0][0]
    bcd101 = prob[0][1][0][1] + prob[1][1][0][1]
    bcd110 = prob[0][1][1][0] + prob[1][1][1][0]
    bcd111 = prob[0][1][1][1] + prob[1][1][1][1]

    # Compute the sum of x * log(x) for each term
    result = (abc000 * np.log(abc000) + abc001 * np.log(abc001) + abc010 * np.log(abc010) + abc011 * np.log(abc011) +
    abc100 * np.log(abc100) + abc101 * np.log(abc101) + abc110 * np.log(abc110) + abc111 * np.log(abc111) +
    abd000 * np.log(abd000) + abd001 * np.log(abd001) + abd010 * np.log(abd010) + abd011 * np.log(abd011) +
    abd100 * np.log(abd100) + abd101 * np.log(abd101) + abd110 * np.log(abd110) + abd111 * np.log(abd111) +
    acd000 * np.log(acd000) + acd001 * np.log(acd001) + acd010 * np.log(acd010) + acd011 * np.log(acd011) +
    acd100 * np.log(acd100) + acd101 * np.log(acd101) + acd110 * np.log(acd110) + acd111 * np.log(acd111) +
    bcd000 * np.log(bcd000) + bcd001 * np.log(bcd001) + bcd010 * np.log(bcd010) + bcd011 * np.log(bcd011) +
    bcd100 * np.log(bcd100) + bcd101 * np.log(bcd101) + bcd110 * np.log(bcd110) + bcd111 * np.log(bcd111))

    return result*para

def compute_point_cluster_energy(prob,lattice="BCC"):
    if lattice=="BCC":
        para=-0.25
    elif lattice=="FCC":
        para=1.25
    a0 = prob[0][0][0][0] + prob[0][0][0][1] + prob[0][0][1][0] + prob[0][0][1][1] + prob[0][1][0][0] + prob[0][1][0][1] + prob[0][1][1][0] + prob[0][1][1][1]
    a1 = prob[1][0][0][0] + prob[1][0][0][1] + prob[1][0][1][0] + prob[1][0][1][1] + prob[1][1][0][0] + prob[1][1][0][1] + prob[1][1][1][0] + prob[1][1][1][1]
    b0 = prob[0][0][0][0] + prob[0][0][0][1] + prob[0][0][1][0] + prob[0][0][1][1] + prob[1][0][0][0] + prob[1][0][0][1] + prob[1][0][1][0] + prob[1][0][1][1]
    b1 = prob[0][1][0][0] + prob[0][1][0][1] + prob[0][1][1][0] + prob[0][1][1][1] + prob[1][1][0][0] + prob[1][1][0][1] + prob[1][1][1][0] + prob[1][1][1][1]
    c0 = prob[0][0][0][0] + prob[0][0][0][1] + prob[0][1][0][0] + prob[0][1][0][1] + prob[1][0][0][0] + prob[1][0][0][1] + prob[1][1][0][0] + prob[1][1][0][1]
    c1 = prob[0][0][1][0] + prob[0][0][1][1] + prob[0][1][1][0] + prob[0][1][1][1] + prob[1][0][1][0] + prob[1][0][1][1] + prob[1][1][1][0] + prob[1][1][1][1]
    d0 = prob[0][0][0][0] + prob[0][0][1][0] + prob[0][1][0][0] + prob[0][1][1][0] + prob[1][0][0][0] + prob[1][0][1][0] + prob[1][1][0][0] + prob[1][1][1][0]
    d1 = prob[0][0][0][1] + prob[0][0][1][1] + prob[0][1][0][1] + prob[0][1][1][1] + prob[1][0][0][1] + prob[1][0][1][1] + prob[1][1][0][1] + prob[1][1][1][1]
    result = (a0 * np.log(a0) + a1 * np.log(a1) +b0 * np.log(b0) + b1 * np.log(b1) +
    c0 * np.log(c0) + c1 * np.log(c1) +d0 * np.log(d0) + d1 * np.log(d1))
    point_prob=np.array([[a0,b0,c0,d0],[a1,b1,c1,d1]])

    return result*para,point_prob

# Compute the sum of x * log(x) for each term
    

class userinput:             #read in user input
    def __init__(self,usrinput):
        inputdictionary={
            "-t":False,
            "-shell":False,
            "main.py":True,
            "-c":[2,"int"],
            "-maxsize":[4,"int"],
            "-cs":["FCC","string"],
            "-basic":["basic.in","string"],
            "-point":False,
            "-vib":[1.0,"float"],
            "-e":[0,"float"],
            "-h":False,
            "-disp":False,
            "-name":["phasediagram.png","string"],
            "-dmu":[0.05,"float"],
            "-unit":["hypo","string"],
            "-control":["control.txt","string"],
            "-e12":[0.0,"float"],
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
        self.vibration_parameter=self.inputdictionary['-vib'][0]
        self.displayhelptext=self.inputdictionary['-h']
        self.elastic_parameter=self.inputdictionary['-e'][0]
        self.phasediagram_name=self.inputdictionary['-name'][0]
        #potential_precision=inputs.inputdictionary['-dmu'][0]
        self.unituse=self.inputdictionary['-unit'][0]
        self.controlfilename=self.inputdictionary['-control'][0]
        self.first_second_bondratio=self.inputdictionary['-e12'][0]
        self.read_cluster_energy()
    
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
    
    def read_cluster_energy(self):
        txt=np.loadtxt(self.inputdictionary['-basic'][0])
        if txt.ndim==1:
            self.basic_cluster_energy,self.local_energy_parameter=txt,np.ones(np.shape(txt)[0])
        else:
            self.basic_cluster_energy,self.local_energy_parameter=txt[0],txt[1]

class phaseboundary:
    def __init__(self,Tstart,mustart,muend,lowphase,highphase,c1,c2,direction=1):
        self.status="open"
        self.signal=0
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

def compare_phb(phb1,phb2,mutol=0.01,Ttol=0.005):   #phb1 is the phb to be determined and phb2 is the existing phb
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

#need editing
def order_phb(phb):
    if phb[2][1]<phb[2][0]:
        return np.flip(phb,axis=1)
    else:
        return phb

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
                "Tmax":5.0,
                "Tstart":2.0,
                "R":1.0,
                "mutol":0.0101,
                "mumin":0.0,
                "mumax":100.0
            }
            return control_dictionary
        elif unit=="bcc":
            control_dictionary={
                "dT":0.2,             
                "dmurough":0.2,
                "dmuprecise":0.02,
                "dmuscan":1.0,
                "Tmax":30.0,
                "Tstart":2.0,
                "R":1.0,
                "mutol":0.0101,
                "mumin":0.0,
                "mumax":100.0
            }
            return control_dictionary        
        elif unit=="kj":
            control_dictionary={              #adjust later for real unit system
                "dT":5.0,             
                "dmurough":0.1,
                "dmuprecise":0.01,
                "dmuscan":0.2,
                "Tmax":520.0,
                "Tstart":490.0,           #change to 300 later
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
    #now it is BCC
    elif phase=="B2":
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

def edit_local_parameter(index,value,file_path="basic.in"):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    second_row = lines[1].strip().split()
    second_row[index] = str(round(value,2))  # Replace 'new_value' with the value you want
    # Join the elements back into a string
    lines[1] = ' '.join(second_row) + '\n'
    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

def write_invariant_point(coordinate, logfilename="log.txt"):
    with open(logfilename, 'a') as logfile:
        logfile.write(f"Invariantpoint: {coordinate[0]} {coordinate[1]}\n")

def write_eutectic_point(coordinate, logfilename="log.txt"):
    with open(logfilename, 'a') as logfile:
        logfile.write(f"Eutecticpoint: {coordinate[0]} {coordinate[1]}\n")

if __name__ == '__main__':
    coordinate = [1.2, 3.4]
    write_invariant_point(coordinate)
    write_eutectic_point(coordinate)
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