#!/usr/bin/python

import numpy as np
import os
import sys
from fylcvm import FYLCVM
from fylcvmfcc import FCC
import utility
import time
import warnings

def CVM_loader(structure,number_of_component,max_clustersize,basic_cluster_energy,vibration_parameter,local_energy_parameter,elastic_parameter,control_dict):
    #need a function that detects all existing methods
    crystal_structures={
        "fcc":FCC,
        #"bcc":BCC,
        "triclinic":FYLCVM
    }
    myCVM=crystal_structures.get(structure.lower(),FYLCVM)
    return myCVM(number_of_component,max_clustersize,basic_cluster_energy,vibration_parameter,local_energy_parameter,elastic_parameter,control_dict)

if __name__ == '__main__':
    print("Welcome to FYL-CVM code")
    warnings.filterwarnings('ignore')
    start_time = time.time()
    argv=["-t"]
    #argv=["-t","-unit=kj","-vib=1.2"]
    #argv=["-t","-vib=1.1"]
    #inputs=utility.userinput(sys.argv)              
    inputs=utility.userinput(argv)
    inputs.read_input()

    if inputs.displayhelptext==True:
        utility.displayhelptext()
    control_dict=utility.read_control_variable(inputs.controlfilename,inputs.unituse)
    #test the code
    if inputs.testmode:
        task="compute"
        composition=[0.38225,0.61775]
    
        print("code test mod")
        myCVM=CVM_loader("FCC",inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy
                         ,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)
        #myCVM=FYLCVM(inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)                                #maybe considering using a class creator later
        if task=="point":
            print("compute the energy of given point")
            result=myCVM.optimize_free_energy(1,composition,"basinhopping")
            myCVM.output_optimization(result,1,composition)
        elif task=="trace":
            dictin={'range': [7.4,7.6], 'T': 1, 'phase': ['L10', 'L12'],'composition':[0.42,0.38]}
            phb=utility.phaseboundary(1,7.4,7.6,'L10', 'L12',0.42,0.38)
            #(newphb,muT,starting_pointnumber)=myCVM.trace_phase_boundary_v1(dictin)
            newphb=myCVM.trace_phase_boundary_v2(phb)
            print(phb.x1mat)
            print(phb.x2mat)
            print(phb.Tspace)
            #(newphb,muT,starting_pointnumber)=myCVM.trace_phase_boundary(dictin)
            #print(myCVM.starting_point_list)
            print(len(myCVM.node_list))
        elif task=="compute":
            #myCVM.compute_phase_diagram(control_dict['Tstart'],0,25,inputs.phasediagram_name)
            myCVM.compute_phase_diagram_v2(control_dict['Tstart'],0,10,inputs.phasediagram_name)
            #myCVM.compute_phase_diagram(2.6,0,10,inputs.phasediagram_name)
            #print(myCVM.phase_boundary_list)
            #myCVM.plot_phase_diagram0(output_name)
        elif task=="scatter":
            #myCVM.scan_phase_diagram(control_dict['Tstart'],-75,75)
            myCVM.scan_phase_diagram(1,0,10)
        elif task=="scan":
            myCVM.scan_phase_boundary(1,0,10)
        elif task=="find":
            myCVM.find_phase_boundary(500,900,[0.5,0.5],"brute")
        elif task=="brute":
            inputlist=[(10,-20,"L10"),(10,-10,"L10"),(10,0,"L12"),(10,10,"L12")]
            for inputvalue in inputlist:
                myCVM.plot_potential_surface(inputvalue[0],inputvalue[1],inputvalue[2])
        elif task=="test":
            print("just check if input is correct")
            print(inputs.phasediagram_name)
        print("--- %s seconds ---" % (time.time() - start_time))
    elif inputs.shellmode:
        print("shell mode in construction")
