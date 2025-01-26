#!/usr/bin/python

import numpy as np
import os
import sys
from fylcvm import FYLCVM
from fylcvmfcc import FCC
from fylcvmbcc import BCC
import utility
import time
import warnings

def CVM_loader(structure,inputs,control_dict):
    #need a function that detects all existing methods
    crystal_structures={
        "fcc":FCC,
        "bcc":BCC,
        "triclinic":FYLCVM
    }
    myCVM=crystal_structures.get(structure.lower(),FYLCVM)
    return myCVM(inputs,control_dict)

if __name__ == '__main__':
    print("Welcome to FYL-CVM code")
    warnings.filterwarnings('ignore')
    start_time = time.time()
    #argv=["-t","-unit=kj","-basic=real.in"]
    #argv=["-t","-unit=kj","-vib=1.2"]
    argv=["-t","-e12=0.0","-unit=bcc"]
    #inputs=utility.userinput(sys.argv)              
    inputs=utility.userinput(argv)
    inputs.read_input()

    if inputs.displayhelptext:
        utility.displayhelptext()
    control_dict=utility.read_control_variable(inputs.controlfilename,inputs.unituse)
    #test the code
    if inputs.testmode:
        task="compute"
        composition=[0.38225,0.61775]
    
        print("code test mod")
        myCVM=CVM_loader("BCC",inputs,control_dict)

        if task=="compute":
            #myCVM.compute_phase_diagram(control_dict['Tstart'],0,25,inputs.phasediagram_name)
            myCVM.compute_phase_diagram_v2(control_dict['Tstart'],0,100,inputs.phasediagram_name)
        elif task=="scatter":
            #myCVM.scan_phase_diagram(control_dict['Tstart'],-75,75)
            myCVM.scan_phase_diagram(1,0,10)
        elif task=="scan":
            myCVM.scan_phase_boundary(1,0,10)
        elif task=="find":
            myCVM.find_phase_boundary(500,900,[0.5,0.5],"brute")
        elif task=="area":
            #mu is between 24.23 and 24.25 composition is between 0.408 and 0.389 at T=500.00
            Trange=np.array([490,520])
            murange=np.array([23,45])
            myCVM.plot_point_region(Trange,murange,5,0.1)
        elif task=="test":
            print("just check if input is correct")
            print(myCVM.basicclusterenergy)
        print("--- %s seconds ---" % (time.time() - start_time))
    elif inputs.shellmode:
        '''
        print("write whatever you want here")
        myCVM=CVM_loader("FCC",inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy
                         ,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)
        
        T=1.2
        c1mat=np.array([])
        c2mat=np.array([])
        for i in range(16):
            myCVM.update_parameter("e",i)
        #print(myCVM.elastic_parameter)
            mustart,muend,c1,c2,Estart,Eend,result,result2=myCVM.search_phase_boundary_v1(T,3,0.5,steps=20)
            c1mat=np.append(c1mat,min(c1,c2))
            c2mat=np.append(c2mat,max(c1,c2))
        name="elastic_c"+str(T)+".txt"
        f=open(name,"w")
        print([np.vstack((c1mat,c2mat))],file=f)
        f.close()
        utility.replace_word_in_file(name,"array","np.array") '''    
        #work on Jan 2nd
        '''
        print("study the effect of local parameter")
        for i in range(3):
            p=1.05+0.02*i
            name="local"+str(round(p,2))+".txt"
            phasediagram_name="local"+str(round(p,2))+".png"
            utility.edit_local_parameter(2,p)             #update file first
            inputs.read_cluster_energy()                             #then update input class
            myCVM=CVM_loader("FCC",inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy
                         ,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)
            myCVM.compute_phase_diagram_v2(control_dict['Tstart'],0,15,phasediagram_name,name)
            del myCVM'''
        #work on Jan 3rd
        print("study the effect of elastic parameter")
        for i in range(1,15):
            phasediagram_name="elastic"+str(i)+".png"
            name="elastic"+str(i)+".txt"
            myCVM=CVM_loader("FCC",inputs,control_dict)
            myCVM.update_parameter("e",i)
            myCVM.compute_phase_diagram_v2(control_dict['Tstart'],0,10,phasediagram_name,name)
            del myCVM
        
        #myCVM.plotGx(1.0,0,25,0.1)
        print("--- %s seconds ---" % (time.time() - start_time))