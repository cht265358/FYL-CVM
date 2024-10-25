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
    #argv=["-t","-vib=0.95"]
    #argv=["-t","-unit=kj"]
    argv=["-t"]
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
    task="compute"
    composition=[0.38225,0.61775]

    print("code test mod")
    myCVM=FYLCVM(number_of_component,max_clustersize,basic_cluster_energy,vibration_parameter,elastic_parameter,control_dict)                                #maybe considering using a class creator later
    if task=="point":
        print("compute the energy of given point")
        result=myCVM.optimize_free_energy(1,composition,"basinhopping")
        myCVM.output_optimization(result,1,composition)
    elif task=="trace":
        dictin={'range': [18.2,18.4], 'T': 600, 'phase': ['L10', 'A1']}
        (newphb,starting_pointnumber)=myCVM.trace_phase_boundary(dictin)
        print(myCVM.starting_point_list)
    elif task=="compute":
        myCVM.compute_phase_diagram(control_dict['Tstart'],0,30)
        print(myCVM.phase_boundary_list)
        #myCVM.plot_phase_diagram0(output_name)
    elif task=="scan":
        myCVM.scan_phase_boundary(650,0,20)
    elif task=="find":
        myCVM.find_phase_boundary(500,900,[0.5,0.5],"brute")
    elif task=="brute":
        inputlist=[(10,-20,"L10"),(10,-10,"L10"),(10,0,"L12"),(10,10,"L12")]
        for inputvalue in inputlist:
            myCVM.plot_potential_surface(inputvalue[0],inputvalue[1],inputvalue[2])
    
    print("--- %s seconds ---" % (time.time() - start_time))


    

'''
        inputdictionary={
            "-t":[False,"bool"],
            "main.py":[True,"bool"],
            "-c":[2,"int"],
            "-maxsize":[4,"int"],
            "-cs":["FCC","string"]
        }'''

