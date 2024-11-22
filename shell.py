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
    print("Here is the shell workspace for FYL-CVM code")
    warnings.filterwarnings('ignore')
    start_time = time.time()

    #start to write your code here
    argv=["-t","-vib=0.9"]      
    inputs=utility.userinput(argv)
    inputs.read_input()
    control_dict=utility.read_control_variable(inputs.controlfilename,inputs.unituse)

    myCVM=CVM_loader("FCC",inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy
                        ,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)
    
    myCVM.compute_phase_diagram_v2(control_dict['Tstart'],5,10,inputs.phasediagram_name)
