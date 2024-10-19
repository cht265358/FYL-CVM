#!/usr/bin/python

import numpy as np
import os
import sys
from fylcvm import FYLCVM
import utility
import time
import warnings

if __name__ == '__main__':
    print("Hello FYL-CVM code world")
    warnings.filterwarnings('ignore')
    start_time = time.time()
    argv=["-t","-unit=kj"]
    #inputs=utility.userinput(sys.argv)              
    inputs=utility.userinput(argv)
    inputs.read_input()
    #write in inputs
    number_of_component=inputs.inputdictionary['-c'][0]
    max_clustersize=inputs.inputdictionary['-maxsize'][0]
    crystal_system=inputs.inputdictionary['-cs'][0]
    basic_cluster_energy=np.loadtxt(inputs.inputdictionary['-basic'][0])
    vibration_parameter=inputs.inputdictionary['-vib'][0]
    displayhelptext=inputs.inputdictionary['-h']
    elastic_parameter=inputs.inputdictionary['-e'][0]
    output_name=inputs.inputdictionary['-name'][0]
    #potential_precision=inputs.inputdictionary['-dmu'][0]
    unituse=inputs.inputdictionary['-unit'][0]
    if displayhelptext==True:
        utility.displayhelptext()
    controlfilename=inputs.inputdictionary['-control'][0]
    control_dict=utility.read_control_variable(controlfilename,unituse)
    #test the code
    task="scan"
    composition=[0.38225,0.61775]
    if inputs.inputdictionary['-t']:
        print("code test mod")
        myCVM=FYLCVM(number_of_component,max_clustersize,basic_cluster_energy,vibration_parameter,elastic_parameter,control_dict)                                #maybe considering using a class creator later
        if inputs.inputdictionary['-point']==True:
            print("compute the energy of given point")
            result=myCVM.optimize_free_energy(1,composition,"basinhopping")
            myCVM.output_optimization(result,1,composition)
        else:
            if task=="trace":
                dictin={'range': [7.4, 7.5], 'T': 1.0, 'phase': ['L10', 'L12']}
                (newphb,starting_pointnumber)=myCVM.trace_phase_boundary(dictin)
                print(myCVM.startingpointlist)
            elif task=="compute":
                myCVM.compute_phase_diagram(control_dict['Tstart'],-100,100)
                print(myCVM.phase_boundary_list)
                myCVM.plot_phase_diagram0(output_name)
            elif task=="scan":
                myCVM.scan_phase_boundary(700,-100,100)
    
    print("--- %s seconds ---" % (time.time() - start_time))


    

'''
        inputdictionary={
            "-t":[False,"bool"],
            "main.py":[True,"bool"],
            "-c":[2,"int"],
            "-maxsize":[4,"int"],
            "-cs":["FCC","string"]
        }'''

