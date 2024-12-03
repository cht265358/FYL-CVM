#!/usr/bin/python

#Bragg-Williams example template

import numpy as np
import os
import sys
from scipy.optimize import minimize
from scipy.optimize import basinhopping

class BW:       #BWer,local_parameter,elastic_para,control_dict):                #initialize all variables
    def __init__(self,input_parameters):
        #self.parameter=paramter    #read in input parameters
        print("initialize the class")

    def compute_grand_potential(self,mu,T,additional_parameter):    #compute the grand cononical potential
        #total_energy=BWenergy_expression
        grand_potential=total_energy+composition[1]*mu
        return grand_potential
    
    def compute_grand_potential_output(self,mu,T,additional_parameter):    #compute the grand cononical potential
        #basically the same code with the previous function but return additional variables as minimization usually requires only 1 return variable
        return grand_potential,composition,F,E
    
    def optimize_grand_potential(self,T,mu,method="NM"):   #you can adjust the method based on input string
        guess=np.array([2,2,2,-1])        #adjust the guess based on your variable
        if method=="NM":
            initial_guess2=np.array([[0,0,0,0],[5,5,5,0],[5,5,0,5],[5,0,5,5],[0,5,5,5]])
            options={
                'initial_simplex':initial_guess2,
                'fatol':1.0e-8,
                'xatol':1.0e-6
            }
            positive=((-10,10),(-10,10),(-10,10),(-10,10))
            result=minimize(self.compute_grand_potential,guess,method='Nelder-Mead',options=options,args=(T,mu),bounds=positive)
        elif method=="basinhopping":
            minimizer_kwargs={"args":(T,mu)}
            result=basinhopping(self.compute_grand_potential,guess,minimizer_kwargs=minimizer_kwargs)
        return result        #result is a python built-in class, for more information please check scipy.optimize.minimize documentary
    
    #This section is the more complex trace phase boundary code
    '''
    def trace_phase_boundary_v2(self,phb):           #always move in the direction of incresing mu
        lowphase=phb.lowphase
        highphase=phb.highphase
        mustart=phb.mustart
        T=phb.Tstart
        print("######trace phase boundary between "+lowphase+" and "+highphase+" ######")
        print("lowphase is "+lowphase)
        print("highphase is "+highphase)
        print("T is "+str(T))
        print("direction is "+str(phb.direction))
        count=1

        while (self.Tmin<=T<=self.Tmax):
            T=self.dT*count*phb.direction+phb.Tstart
            mustart+=utility.compute_dmu(phb.muspace)             #edit after compute_dmu() has value
            result=self.optimize_grand_potential(T,mustart,"BFGS_v2")               #try to get rid of this step later
            currentphase=self.identify_phase(result.x)
            #print("current phase is "+currentphase)
            if currentphase==lowphase:
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustart,self.dmurough,20)
                if muend==1000:
                    print("phase boundary ends at assigned boundary")
                    phb.status="boundary"
                    return phb
                elif muend==False:
                    #This corresponds to invariant point so start a new node connected to the old node, also the direction is -1
                    print("Find invariant point1!!!!!!!!1!!!!!!!!!!!!!!!!")
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,muendsave,self.dmurough,20)
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                        return phb
                    else:
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-1)
                        phb.status="exit"
                        return phb        
                if self.identify_phase(result2.x)!=highphase:                  #new phase found
                    print("new phase "+self.identify_phase(result2.x)+" detected, follow two new boundary")
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,muend,self.dmurough,10)
                    if muend:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                        phb.status="exit"
                    else:
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,muendsave,self.dmurough,10)
                        if mustart:
                            self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-1)
                            phb.status="exit"
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                    return phb
            elif currentphase==highphase:
                #search_phase_boundary_v1 always from return mu low to high
                print("highphase is "+highphase)
                print("mustart is "+str(mustart))
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustart,-self.dmurough,20)
                if muend==1000:
                    print("phase boundary ends at assigned boundary")
                    phb.status="boundary"
                    return phb
                elif muend==False:
                    #This corresponds to invariant point so start a new node connected to the old node, also the direction is -1
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,20)
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                        return phb
                    else:
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-1)
                        phb.status="exit"
                        return phb   
                if self.identify_phase(result1.x)!=lowphase:
                    print("new phase "+self.identify_phase(result1.x)+" detected, follow two new boundary")
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustart,-self.dmurough,10)
                    if muend:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                        phb.status="exit"
                    else:
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,10)
                        if mustart:
                            self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-1)
                            phb.status="exit"
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                    return phb
            else:   
                mustartuse=mustart       
                print("new phase "+currentphase+" detected, start two new search")
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,self.dmurough,20)    #one forward and one reverse
                print("mustart is "+str(mustart))
                print("muend is "+str(muend))
                if muend and muend!=1000:   #edit later for unlikely boundary condition
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,-self.dmurough,20)
                    if muend and muend!=1000:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                        phb.status="exit"
                    elif muend==1000:
                        phb.status="exit"
                    elif muend==False:
                        print("we are here?")
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,10)
                        if mustart:
                            self.start_new_phb(self.dT*phb.direction,mustart,muend,result1,result2,x1,x2)
                            phb.status="exit"
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                elif muend==False:
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,self.dmurough,10)    #The reverse search must have boundary
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                    else:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2)
                        phb.status="exit"
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartuse,self.dmurough,10)    #The reverse search must have boundary
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                    else:
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-1)
                        phb.status="exit"
                return phb

            #print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2)+" at T="+str(T))
            print("mu is between {:.2f} and {:.2f} composition is between {:.3f} and {:.3f} at T={:.2f}".format(mustart,muend,x1,x2,T))
            #print("sitepot1 is "+str(result1.x)+" sitepot2 is "+str(result2.x))     
            phb.x1mat=np.append(phb.x1mat,x1)
            phb.x2mat=np.append(phb.x2mat,x2)
            phb.Tspace=np.append(phb.Tspace,T)
            phb.muspace=np.append(phb.muspace,0.5*(mustart+muend))
            mustartsave=mustart
            muendsave=muend
            count+=1
        
        phb.status="exit"
        return phb
    
    def search_phase_boundary_v1(self,T,mustart,steplength,steps=20):
        cdata=np.zeros(steps)
        for i in range(steps):
            muuse=mustart+i*steplength
            if muuse>self.mumax or muuse<self.mumin:
                print("boundary encountered")
                return False,1000,False,False,False,False,False,False
            result=self.optimize_grand_potential(T,muuse,"BFGS_v2")  
            #print("current phase is "+str(self.identify_phase(result.x))+" at potential "+str(muuse)+" para is "+str(result.x))
            (potential,composition,F,E)=self.compute_grand_potential_output(result.x,T,muuse)
            cdata[i]=composition[0]
            currentphase=self.identify_phase(result.x)
            if i>0:
                if currentphase!=phasesave:
                    print("phase transition between "+currentphase+" and "+phasesave+"\nmu from "+str(muuse-steplength)+" to "+str(muuse))          
                    return self.search_phase_boundary_bisection(T,muuse-steplength,muuse,phasesave,currentphase)
            phasesave=currentphase
            Esave=E
            resultsave=result
        print("can't find phase boundary")
        return False,False,cdata[i-1],cdata[i],Esave,E,resultsave,result
    
    def search_phase_boundary_bisection(self,T,mustart,muend,phasestart,phaseend,composition=np.array([0,1])):        #only use this function to search within a given range
        if mustart>muend:
            mustart,muend=muend,mustart         #make sure mustart is always smaller
            phasestart,phaseend=phaseend,phasestart
        resultlist=[]
        while np.abs(mustart-muend)>self.mutol:
            mumid=0.5*(mustart+muend)
            result=self.optimize_grand_potential(T,mumid,"BFGS_v2") 
            resultlist.append(result)
            phasemid=self.identify_phase(result.x)
            if phasemid!=phasestart and phasemid!=phaseend:        #need recursion
                print("new phase "+phasemid+" detected?")
            elif phasemid==phaseend:
                muend=mumid
            elif phasemid==phasestart:
                mustart=mumid
        
        #use logic to remove last run later
        if mustart==mumid:
            (potential1,cstart,F1,Estart)=self.compute_grand_potential_output(result.x,T,mustart)
            result2=self.optimize_grand_potential(T,muend,"BFGS_v2")    #optimize here later
            (potential2,cend,F2,Eend)=self.compute_grand_potential_output(result2.x,T,muend)
        else:
            (potential1,cend,F1,Eend)=self.compute_grand_potential_output(result.x,T,muend)
            result2=self.optimize_grand_potential(T,mustart,"BFGS_v2")    #optimize here later
            (potential2,cstart,F2,Estart)=self.compute_grand_potential_output(result2.x,T,mustart)
            result,result2=result2,result   #swap to make sure mustart match result and muend match result2
        return mustart,muend,cstart[0],cend[0],Estart,Eend,result,result2
        '''

if __name__ == '__main__':
    yourBW=BW(input_parameters)
    BW.optimize_grand_potential(T,mu)