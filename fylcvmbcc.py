#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import plotter
import matplotlib.pyplot as plt
import freeenergy
from scipy.optimize import minimize
from scipy.optimize import brute
from scipy.optimize import basinhopping
from multiprocessing import Pool
from matplotlib import cm
from fylcvm import FYLCVM
#svm stands for site_variable_matrix

def unequal(i,j):
    return 1 if i != j else 0

class BCC(FYLCVM):       #sub class for BCC
    def __init__(self,inputs,control_dict):                #initialize all variables
        super().__init__(inputs,control_dict)
        self.lattice="BCC"
        #The current code only read from parameter (easier)
        self.map_basic_cluster_energy(inputs.first_second_bondratio)

    def map_basic_cluster_energy(self,E):     #E can be either filename or parametervalue
        self.basicclusterenergy=np.zeros(self.shape)                  
        if os.path.exists(str(E)):
            Edata = np.loadtxt(str(E))
            #edit later for TO approximation    
            for i,j,k,l in itertools.product(range(self.component), repeat=4):       
                self.basicclusterenergy[i][j][k][l]=E[i+k+j+l]
        else:
            E1=-1
            E2=E*E1
            for i,j,k,l in itertools.product(range(self.component), repeat=4):#for BCC the 4 positions are distinct 
                self.basicclusterenergy[i][j][k][l]=E2*(unequal(i,j)+unequal(k,l))+E1*(unequal(i,k)+unequal(i,l)+unequal(j,k)+unequal(j,l))
    
    def map_vibrational_energy(self,T):        #Assume symmetry wAA=wBB
        #Fmat=[0,1.5*T*self.R*np.log(self.vib_para*self.local_para[1]),2*T*self.R*np.log(self.vib_para*self.local_para[2])
        #      ,1.5*T*self.R*np.log(self.vib_para*self.local_para[3]),0]
        #Fmat=[0,1.5*T*self.R*np.log(self.vib_para),2*T*self.R*np.log(self.vib_para*1.05),1.5*T*self.R*np.log(self.vib_para),0]
        for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
            #self.vibrationalenergy[i][j][k][l]=Fmat[i+k+j+l]
            self.vibrationalenergy[i][j][k][l]=0

    def compute_partition_function(self,svm,T,Fvib):
        partition=0                                #now skip site variable
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            partition+=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))
        return partition
        
    def compute_basic_cluster_prob(self,partition,svm,T,Fvib):        #svm is site potential matrix
        basic_cluster_prob=np.zeros(self.shape)
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            basic_cluster_prob[i][j][k][l]=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))/partition
        return basic_cluster_prob

    def identify_phase(self,site_potential_input):                           #Figure out D03, B2 and A2 first
        error=1.0e-3
        x1=site_potential_input[0]
        x2=site_potential_input[1]
        x3=site_potential_input[2]
        x4=site_potential_input[3]
        x=0.25*(x1+x2-x3-x4)
        y=0.5*(x3-x4)
        z=0.5*(x1-x2)
        if np.abs(z)<error and np.abs(y)>error:
            return "D03"
        elif np.abs(z)<error and np.abs(y)<=error and np.abs(x)>error:
            return "B2"
        elif np.abs(z)<error and np.abs(y)<=error and np.abs(x)<=error:
            return "A2"

    def compute_grand_potential(self,site_potential_input,T,mu):    #compute the grand cononical potential
        self.map_vibrational_energy(T)
        spm=np.zeros((self.component,self.clustersize))             #spm is site potential matrix
        for i in range(self.clustersize):
            spm[1][i]=site_potential_input[i]
        partition=self.compute_partition_function(spm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,spm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0
        
        three_body_energy=freeenergy.compute_ternary_cluster_energy(prob,"BCC")*self.R*T
        two_body_energy=freeenergy.compute_binary_cluster_energy(prob,"BCC")*self.R*T
        pointenergy,point_prob=freeenergy.compute_point_cluster_energy(prob,"BCC")
        pointenergy=pointenergy*self.R*T

        composition=np.zeros(self.component)
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                composition[j]+=point_prob[j][i]*0.25
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=self.R*T*np.log(partition)
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        
        total_energy=6*basic_cluster_energy+three_body_energy+two_body_energy+pointenergy-composition[1]*mu+elastic_energy
        #entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        #print("energy is "+str(total_energy)+" and para is "+str(site_potential_input))
        #print("basic energy is "+str(6*basic_cluster_energy)+" point energy "+str(pointenergy)+" binary energy "+str(two_body_energy)+"ternary energy "+str(three_body_energy))
        #print("composition is "+str(composition[0]))
        return total_energy
    
    def compute_grand_potential_output(self,site_potential_input,T,mu):    #compute the grand cononical potential
        self.map_vibrational_energy(T)
        spm=np.zeros((self.component,self.clustersize))             #spm is site potential matrix
        for i in range(self.clustersize):
            spm[1][i]=site_potential_input[i]
        partition=self.compute_partition_function(spm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,spm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0

        three_body_energy=freeenergy.compute_binary_cluster_energy(prob,"BCC")*self.R*T
        two_body_energy=freeenergy.compute_binary_cluster_energy(prob,"BCC")*self.R*T
        pointenergy,point_prob=freeenergy.compute_point_cluster_energy(prob,"BCC")
        pointenergy=pointenergy*self.R*T

        composition=np.zeros(self.component)
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                composition[j]+=point_prob[j][i]*0.25
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=self.R*T*np.log(partition)
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        
        total_energy=6*basic_cluster_energy+three_body_energy+two_body_energy+pointenergy-composition[1]*mu+elastic_energy
        F=total_energy-composition[1]*mu
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        return total_energy,composition,F,0
    
    def optimize_grand_potential(self,T,mu,method):
        if method=="BFGS_v2":
            guess_1=np.array([2,-1])
            guess_2=np.array([-1,2])
            l10_1=np.array([0,0,])   
            l10_2=np.array([-5,5,])
            guesslist=(guess_1,guess_2,l10_1,l10_2)  
            result,phase=self.multiBFGS_v2(guesslist,(T,mu))
            result.x=utility.map_variable_to_phase(result.x,phase,self.clustersize)         #change later if everything works perfectly, return result.x with size 4
            #print("at potential mu="+str(mu)+" best phase is "+self.identify_phase(result.x))
        else:
            sys.exit("Error: Optimization method can't be found")
        return result
    
    def multiBFGS_v2(self,guesslist,args):             #input two tuple
        Fmin=10
        argsB2=args+("B2",)
        argslist=(argsB2,argsB2,argsB2,argsB2)
        #print(argslist)
        #longguesslist=guesslist*2
        #print(longguesslist)
        inputlist=zip(guesslist,argslist)               #start with A2 B2 first

        with Pool(processes=8) as pool:              #remember to add this back later!!
            resultlist=pool.starmap(self.minimizerv2,inputlist)
            pool.close()

        for result in resultlist:
            #print("phase is "+result[1])
            if result[0].fun<Fmin:
                Fmin=result[0].fun
                bestresult=result[0]
                phase=result[1]
        
        return bestresult,phase

    def minimizerv2(self,initialvalue,args,method='BFGS'):      #return phase
        result=minimize(self.potential_computerv2,initialvalue,method=method,args=args)
        return result,args[2]

    def potential_computerv2(self,initialvalue,T,mu,phase):
        sitepotential_input=np.zeros(self.clustersize)
        if phase=="B2":
            sitepotential_input[0]=sitepotential_input[1]=initialvalue[0]
            sitepotential_input[2]=sitepotential_input[3]=initialvalue[1]
        elif phase=="D03":
            sitepotential_input[0]=sitepotential_input[1]=initialvalue[0]
            sitepotential_input[2]=initialvalue[1]
            sitepotential_input[3]=initialvalue[2]
        else:
            sys.exit("Phase can't be identified")
        return self.compute_grand_potential(sitepotential_input,T,mu)
    
    def find_phase_boundary(self,Tstart,Tstop,composition,method="BFGS"):
        T=Tstart
        count=0
        while T<Tstop:
            result=self.optimize_free_energy(T,composition,method)
            self.output_optimization(result,T,composition,method)
            T+=self.dT
    
    #trace_phase_boundary returns the phase boundary between two specific phases
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
            print("current phase is "+currentphase)
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
            print("current phase is "+str(self.identify_phase(result.x))+" at potential "+str(muuse)+" para is "+str(result.x)+" T is "+str(T))
            (potential,composition,F,E)=self.compute_grand_potential_output(result.x,T,muuse)
            cdata[i]=composition[0]
            print("current composition is "+str(composition[0]))
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

    def scan_phase_boundary_v2(self,T,mustart,muend):        #scan the phase boundary across a wide 
        if mustart>muend:
            mustart,muend=muend,mustart
        print("scan for phase boundary starting point between mu={:.2f} and {:.2f} at T={:.2f}".format(mustart,muend,T))
        muuse=mustart                          #now only search once
        (mu1,muuse,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,muuse,self.dmuscan,int((muend-muuse)/self.dmuscan)+1)
        if mu1:
            self.start_new_phb(T,mu1,muuse,result1,result2,x1,x2)
        else:
            print("can't find phase boundary, please adjust search range")
        print("scan phase boundary finished")          #add scan algorithm later
    
    def scan_phase_diagram(self,Tstart,mustart,muend,phase_diagram_name="scatter.png"):
        #plot the phase diagram by scanning through the whole space
        print("scan through the whole mu-T space")
        phblist=np.empty((0,2))
        count=0
        T=Tstart
        while T<self.Tmax:
            T=Tstart+count*self.dT
            phb,mustart,muend=self.scan_phase_boundary(T,mustart,muend,False)
            if not mustart:
                print("no more phase boundaries, end searching in space")
                break
            phblist=np.append(phblist,phb,axis=0)
            mustart-=50*self.dmuprecise
            muend+=50*self.dmuprecise
            count+=1
        plotter.plot_scatter_rough(phblist,phase_diagram_name)

    def compute_phase_diagram_v2(self,Tstart,mustart,muend,phase_diagram_name="phasediagram.png",phbfile_name="phase_boundary_list.txt",color='r'):
        #first scan phase boundary starting point
        placeholder=self.scan_phase_boundary_v2(Tstart,mustart,muend)
        print("scan finished, now start search")
        unfinishednode=1
        while unfinishednode:  
            unfinishednode=0
            for i in range(len(self.node_list)):
                if self.node_list[i].status=="open":
                    for j in range(len(self.node_list)):
                        if utility.compare_phb(self.node_list[i],self.node_list[j]):
                            self.node_list[i].status="exit"
                    if self.node_list[i].status=="exit":
                        print("node overlap, move to next node")
                        continue
                    unfinishednode=1
                    self.node_list[i]=self.trace_phase_boundary_v2(self.node_list[i])
                    phblist=np.stack((self.node_list[i].x1mat,self.node_list[i].x2mat,self.node_list[i].Tspace))
                    self.phase_boundary_list.append(utility.order_phb(phblist))
        f=open(phbfile_name,"w")
        print(self.phase_boundary_list,file=f)
        f.close()
        utility.replace_word_in_file("phase_boundary_list.txt","array","np.array")          #make the output file that can be directly used by plotter
        plotter.plot_phase_diagram_rough(self.phase_boundary_list)

    def plot_phase_diagram0(self,filename):
        for i in range(len(self.phase_boundary_list)):
            listuse=self.phase_boundary_list[i]
            plt.plot(listuse[0],listuse[2])
            plt.plot(listuse[1],listuse[2])
        plt.savefig(filename)
    
    def start_new_phb(self,Tstart,mustart,muend,result1,result2,x1,x2,direction=1):
        self.node_list.append(utility.phaseboundary(Tstart,mustart,muend,self.identify_phase(result1.x),self.identify_phase(result2.x),x1,x2,direction))

    def output_optimization(self,result,T,component_comp,method):
        if method=="brute":
            result=np.append(result,result[0])
            svm=self.compute_site_variable_matrix(result,component_comp,T) 
            print(" site variable is "+str(svm[1])+" composition is "+str(component_comp)+" free energy is "+str(self.compute_total_energy_brute(result,T,component_comp))+" at T="+str(T))
        else:
            svm=self.compute_site_variable_matrix(result.x,component_comp,T) 
            print("free energy is "+str(result.fun)+" site variable is "+str(svm[1])+" composition is "+str(component_comp))
        print("the phase is "+self.identify_phase(svm[1]))
        return self.identify_phase(svm[1])

    def output_potential(self,result,T,component_comp):
        svm=np.zeros(self.clustersize)
        for i in range(self.clustersize):
            svm[i]=np.exp(result.x[i]/T)
        #print("free energy is "+str(result.fun)+" site variable is "+str(svm[1])+" at temperature "+str(T))
        print("free energy is "+str(result.fun+mu*component_comp[1])+" site variable is "+str(svm)+" composition is "+str(component_comp)+" potential is "+str(mu))

    def print_output(self):
        print("hello world")

if __name__ == '__main__':
    print("you are in FYLCVM code?")