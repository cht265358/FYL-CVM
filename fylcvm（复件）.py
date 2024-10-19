#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import brute
from scipy.optimize import basinhopping
from multiprocessing import Pool
#svm stands for site_variable_matrix

class FYLCVM:       #base class for FCC
    def __init__(self,component,maxsize,E,vibration_parameter,elastic_para,potential_precision,Rin):                #initialize all variables
        print("start FYL-CVM")
        self.Ruse=Rin
        print("R is "+str(self.Ruse))
        N=6.02e23;       #number of particles, the size of the system, here we use 1 mole
        self.lattice="FCC"
        self.component=component
        self.clustersize=maxsize
        self.shape=[self.component]*self.clustersize    #ignore using basic cluster with different size first
        self.map_basic_cluster_energy(E)
        self.vibrationalenergy=np.zeros(self.shape)
        self.vibration_parameter=vibration_parameter             #vibration parameter with default of 1
        self.elastic_parameter=elastic_para
        self.starting_point_list=[]
        self.phase_boundary_list=[]
        self.potential_precision=potential_precision

    def map_basic_cluster_energy(self,E):     #map the basic cluster energy to high dimensional array
        self.basicclusterenergy=np.zeros(self.shape)                  #use high dimensional array, input by hand now
        if self.clustersize==4:
            for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
                self.basicclusterenergy[i][j][k][l]=E[i+k+j+l]
    
    '''def map_vibrational_energy(self,wmat,T):        #w is w matrix normalized
        Fmat=[-3*T*np.log(T/wmat[0]),-1.5*T*np.log(T/wmat[0])-1.5*T*np.log(T/wmat[1]),-0.5*T*np.log(T/wmat[0])
             -2*T*np.log(T/wmat[1])-0.5*T*np.log(T/wmat[2]),-1.5*T*np.log(T/wmat[1])-1.5*T*np.log(T/wmat[2]),-3*T*np.log(T/wmat[2])]
        for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
            self.vibrationalenergy[i][j][k][l]=Fmat[i+k+j+l]'''
    
    def map_vibrational_energy(self,vibpara,T):        #Assume symmetry wAA=wBB
        Fmat=[0,-0.75*self.Ruse*T*np.log(vibpara),-self.Ruse*T*np.log(vibpara),-0.75*self.Ruse*T*np.log(vibpara),0]
        for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
            self.vibrationalenergy[i][j][k][l]=Fmat[i+k+j+l]

    def compute_site_variable_matrix(self,site_variable_input,composition,T):         #extend it to larger system if expression is known
        svm=np.ones((self.component,self.clustersize))
        for i in range(self.clustersize-1):
            svm[1][i]=np.exp(site_variable_input[i])                                  #the input is the log of site variable
        nominator=0
        denominator=0
        shape1=[self.component]*(self.clustersize-1)
        for i,j,k in itertools.product(range(self.component), repeat=3):
            nominator+=(composition[1]-(1/self.clustersize)*(i+j+k+0))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][0]/T)
            denominator+=(-composition[1]+(1/self.clustersize)*(i+j+k+1))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][1]/T)
        svm[1][3]=nominator/denominator
        return svm                                                        #check later to adjust for real unit
        
        #merge vibrational energy into basic cluster energy
    def compute_partition_function(self,svm,T,Fvib):
        if svm[0][0]==1:                 #matrix is for site variable
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp((-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/T)  #no need to adjust here
            return partition
        else:                            #matrix is for site potential, svm is actually spm
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/T)
            return partition
        
    def compute_basic_cluster_prob(self,partition,svm,T,Fvib):        #svm is site potential matrix
        basic_cluster_prob=np.zeros(self.shape)
        if svm[0][0]==1:                 #matrix is for site variable
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                basic_cluster_prob[i][j][k][l]=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp((-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.Ruse*T))/partition
        else:                            #matrix is for site potential, svm is actually spm
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                basic_cluster_prob[i][j][k][l]=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/T)/partition
        return basic_cluster_prob

    def identify_phase(self,site_potential_input):                           #identify which phase is showing up
        count=0
        for i in range(self.clustersize):
            if site_potential_input[i]==0:
                return "A1"
            for j in range(i+1,self.clustersize):
                if np.abs((site_potential_input[i]-site_potential_input[j])/site_potential_input[i])<1.0e-2:
                    count+=1
        if count==6:
            return "A1"
        elif count==3:
            return "L12"
        elif count==2:
            return "L10"
        elif count==1:
            return "L12+L10"
            
    def compute_total_energy(self,site_variable_input,T,component_comp):    #site_variable_input is 1d for binary for now
        svm=self.compute_site_variable_matrix(site_variable_input,component_comp,T) 
        partition=self.compute_partition_function(svm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,svm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0
        #two body cluster
        two_body_energy=0
        for i in range(self.clustersize):                 #4*3*2*2 two body cluster
            for j in range(i+1,self.clustersize):
                for k,l in itertools.product(range(self.component), repeat=2):     #decorate this line later
                    position_type_matrix=np.zeros((2,2))
                    position_type_matrix[0][0]=i
                    position_type_matrix[0][1]=j
                    position_type_matrix[1][0]=k
                    position_type_matrix[1][1]=l
                    probij=utility.get_cluster_prob(prob,position_type_matrix)
                    two_body_energy+=probij*np.log(probij)*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=T*np.log(partition)
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        return total_energy
    
    def compute_total_energy_output(self,site_variable_input,T,component_comp):    #simply repeating the function but allow for more output
        svm=self.compute_site_variable_matrix(site_variable_input,component_comp,T) 
        partition=self.compute_partition_function(svm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,svm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0
        #two body cluster
        two_body_energy=0
        for i in range(self.clustersize):                 #4*3*2*2 two body cluster
            for j in range(i+1,self.clustersize):
                for k,l in itertools.product(range(self.component), repeat=2):     #decorate this line later
                    position_type_matrix=np.zeros((2,2))
                    position_type_matrix[0][0]=i
                    position_type_matrix[0][1]=j
                    position_type_matrix[1][0]=k
                    position_type_matrix[1][1]=l
                    probij=utility.get_cluster_prob(prob,position_type_matrix)
                    two_body_energy+=probij*np.log(probij)*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=T*np.log(partition)
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        mu=0
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            product=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]
            mu+=product*np.log(product)*np.exp(self.basicclusterenergy[i][j][k][l])*T

        return total_energy,entropy,mu/partition

    def compute_grand_potential(self,site_potential_input,T,mu):    #compute the grand cononical potential
        #self.map_vibrational_energy(T)
        spm=np.zeros((self.component,self.clustersize))             #spm is site potential matrix
        for i in range(self.clustersize):
            spm[1][i]=site_potential_input[i]
        partition=self.compute_partition_function(spm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,spm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0

        two_body_energy=0.0
        for i in range(self.clustersize):                 #4*3*2*2 two body cluster
            for j in range(i+1,self.clustersize):
                for k,l in itertools.product(range(self.component), repeat=2):     #decorate this line later
                    position_type_matrix=np.zeros((2,2))
                    position_type_matrix[0][0]=i
                    position_type_matrix[0][1]=j
                    position_type_matrix[1][0]=k
                    position_type_matrix[1][1]=l
                    probij=utility.get_cluster_prob(prob,position_type_matrix)
                    two_body_energy+=probij*np.log(probij)*T
        #point cluster  
        pointenergy=0.0
        point_prob=np.zeros((self.component,self.clustersize))
        composition=np.zeros(self.component)
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T
                composition[j]+=point_prob[j][i]*0.25
        #basic cluster start with tetra as default        
        basic_cluster_energy=0.0                          
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=T*np.log(partition)
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy-composition[1]*mu+elastic_energy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T     ##check here
        return total_energy
    
    def compute_grand_potential_output(self,site_potential_input,T,mu):    #compute the grand cononical potential
        #self.map_vibrational_energy(T)
        spm=np.zeros((self.component,self.clustersize))             #spm is site potential matrix
        for i in range(self.clustersize):
            spm[1][i]=site_potential_input[i]
        partition=self.compute_partition_function(spm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,spm,T,self.vibrationalenergy)

        if partition<0:
            return 1000.0

        two_body_energy=0
        for i in range(self.clustersize):                 #4*3*2*2 two body cluster
            for j in range(i+1,self.clustersize):
                for k,l in itertools.product(range(self.component), repeat=2):     #decorate this line later
                    position_type_matrix=np.zeros((2,2))
                    position_type_matrix[0][0]=i
                    position_type_matrix[0][1]=j
                    position_type_matrix[1][0]=k
                    position_type_matrix[1][1]=l
                    probij=utility.get_cluster_prob(prob,position_type_matrix)
                    two_body_energy+=probij*np.log(probij)*T
        #point cluster  
        pointenergy=0
        xb=0
        point_prob=np.zeros((self.component,self.clustersize))
        composition=np.zeros(self.component)
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T
                composition[j]+=point_prob[j][i]*0.25
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=T*np.log(partition)

        E=0
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            E+=self.basicclusterenergy[i][j][k][l]*prob[i][j][k][l]
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy-composition[1]*mu
        F=total_energy-composition[1]*mu
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T  ##check here
        return total_energy,composition,F,E
    
    def optimize_grand_potential(self,T,mu,method):
        guess=np.array([2,2,2,-1])
        if method=="NM":
            initial_guess2=np.array([[0,0,0,0],[5,5,5,0],[5,5,0,5],[5,0,5,5],[0,5,5,5]])
            options={
                'initial_simplex':initial_guess2,
                'fatol':1.0e-8,
                'xatol':1.0e-6
            }
            positive=((-10,10),(-10,10),(-10,10),(-10,10))
            #result=minimize(CVM.compute_total_energy,guess,method='Nelder-Mead',args=(T,component_comp),bounds=positive,tol=1.0e-6)
            result=minimize(self.compute_grand_potential,guess,method='Nelder-Mead',options=options,args=(T,mu),bounds=positive)
        elif method=="basinhopping":
            minimizer_kwargs={"args":(T,mu)}
            result=basinhopping(self.compute_grand_potential,guess,minimizer_kwargs=minimizer_kwargs)
        elif method=="BFGS":         #default value now
            A1=np.array([0,0,0,0])
            L12=np.array([2.25,2.25,2.25,-1.6])
            L10=np.array([3.5,3.5,-0.25,-0.25])     
            guesslist=(A1,L12,L10)
            #result=minimize(self.compute_grand_potential,guess,method='BFGS',args=(T,mu))
            result=self.multiBFGS(guesslist,(T,mu))
        return result
    
    def multiBFGS(self,guesslist,args):             #input two tuple
        Fmin=1000
        argslist=[args]*len(guesslist)
        inputlist=zip(guesslist,argslist)
        with Pool(processes=4) as pool:
            resultlist=pool.starmap(self.minimizer,inputlist)
            pool.close()
            #pool.join()

        for result in resultlist:
            if result.fun<Fmin:
                Fmin=result.fun
                bestresult=result
        return bestresult
    
    def minimizer(self,initialvalue,args,method='BFGS'):      #use a minimizer that can be called by pool.starmap
        result=minimize(self.compute_grand_potential,initialvalue,method=method,args=args)
        return result
    
    '''
    def multiBFGS(self,guesslist,args):             #input two tuple
        Fmin=1000
        for guess in guesslist:
            result=minimize(self.compute_grand_potential,guess,method='BFGS',args=args)
            if result.fun<Fmin:
                Fmin=result.fun
                bestresult=result
        return bestresult'''

    def optimize_free_energy(self,T,component_comp,method):         #main optimization function, optimize site variables at given T and composition
        guess=np.array([2.26,2.26,-1.84])
        #initial_guess=np.ones((self.clustersize-1)*(self.component-1))
        if method=="NM":
            initial_guess2=np.array([[3,3,-1.5],[2.5,2.5,-2],[2.5,2.5,-1.5],[2.64,2.64,-1.897]])
            options={
                'initial_simplex':initial_guess2
            }
            positive=((-20,20),(-20,20),(-20,20))
            #result=minimize(CVM.compute_total_energy,guess,method='Nelder-Mead',args=(T,component_comp),bounds=positive,tol=1.0e-6)
            result=minimize(self.compute_total_energy,guess,method='Nelder-Mead',options=options,args=(T,component_comp),bounds=positive,tol=1.0e-6)
        elif method=="basinhopping":
            minimizer_kwargs={"args":(T,component_comp)}
            result=basinhopping(self.compute_total_energy,guess,minimizer_kwargs=minimizer_kwargs)
        elif method=="brute":
            boundary=(slice(1.5,2.5,0.05),slice(1.5,2.5,0.05),slice(1.5,2.5,0.05))
            result=brute(self.compute_total_energy,boundary,args=(T,component_comp),workers=-1,full_output=False)

        return result
    
    #trace_phase_boundary returns the phase boundary between two specific phases
    def trace_phase_boundary(self,startingdict,steplength=0.05,Tmax=2.5):           #always move in the direction of incresing mu
        lowphase=startingdict['phase'][0]
        highphase=startingdict['phase'][1]
        Tmin=startingdict['T']
        mustart=startingdict['range'][0]
        print("trace phase boundary between "+lowphase+" and "+highphase)
        x1mat=np.array([])
        x2mat=np.array([])
        Tspace=np.array([])               #flexable value
        T=Tmin
        count=0

        while (T<=Tmax):
            T=steplength*count+Tmin
            mustart+=self.compute_dmu()             #edit after compute_dmu() has value
            result=self.optimize_grand_potential(T,mustart,"BFGS")
            currentphase=self.identify_phase(result.x)
            print("current phase is "+currentphase)
            if currentphase==lowphase:
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,self.potential_precision,100)
                if muend==1000:
                    print("phase boundary ends")
                    return np.stack((x1mat,x2mat,Tspace)),0
                if self.identify_phase(result2.x)!=highphase:
                    print("new phase detected, follow two new boundary")
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(newstartingdict)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,0.1,100)
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(newstartingdict)
                    return np.stack((x1mat,x2mat,Tspace)),2
            elif currentphase==highphase:
                #search in reverse so that result is also in reverse
                (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-self.potential_precision,100)
                if muend==1000:
                    print("phase boundary ends")
                    return np.stack((x1mat,x2mat,Tspace)),0
                if self.identify_phase(result1.x)!=lowphase:
                    print("new phase detected, follow two new boundary")
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(newstartingdict)
                    (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-0.1,100)
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(newstartingdict)
                    return np.stack((x1mat,x2mat,Tspace)),2
            else:
                print("new phase detected, start two new search")
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,0.1,100)    #one forward and one reverse
                newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                self.starting_point_list.append(newstartingdict)
                (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-0.1,100)
                newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                self.starting_point_list.append(newstartingdict)
                return np.stack((x1mat,x2mat,Tspace)),2
            #search with higher precision, we can add a loop later
            #(mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart-0.001,0.001,100)
            print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2)+" at T="+str(T))
            print("sitepot1 is "+str(result1.x)+" sitepot2 is "+str(result2.x))     
            x1mat=np.append(x1mat,x1)
            x2mat=np.append(x2mat,x2)
            Tspace=np.append(Tspace,T)
            count+=1
        
        return np.append(np.append(x1mat,x2mat,axis=0),Tspace,axis=0),1
        #plt.plot(x1mat,Tspace)
        #plt.plot(x2mat,Tspace)
        #plt.show()   

    def compute_dmu(self):
        return 0

    def search_phase_boundary(self,T,mustart,steplength,steps):
        if steps>1000:            #max of 1000 steps
            steps=1000
        cdata=np.zeros(steps)
        for i in range(steps):
            muuse=mustart+i*steplength
            #if muuse<0:
                #print("can't find phase boundary")
                #return 1000,1000,cdata[i-1],cdata[i],Esave,E,resultsave,result
            result=self.optimize_grand_potential(T,muuse,"BFGS")
            (potential,composition,F,E)=self.compute_grand_potential_output(result.x,T,muuse)
            cdata[i]=composition[0]
            currentphase=self.identify_phase(result.x)
            if i>0:
                if currentphase!=phasesave:
                    print("phase transition between "+currentphase+" and "+phasesave+"\nmu from "+str(muuse-steplength)+" to "+str(muuse)) 
                    slope=(muuse-0.5*steplength)/T+(E-Esave)/((cdata[i]-cdata[i-1])*T)
                    print("slope is "+str(slope))             
                    return muuse-steplength,muuse,cdata[i-1],cdata[i],Esave,E,resultsave,result
            phasesave=currentphase
            Esave=E
            resultsave=result
        print("can't find phase boundary")
        return 1000,1000,cdata[i-1],cdata[i],Esave,E,resultsave,result
    
    def scan_phase_boundary(self,T,mustart,muend,steplength):        #scan the phase boundary across a wide 
        muuse=mustart
        while muuse<muend:                             #continue search with new muuse until boundary is reached
            (mu1,muuse,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,muuse,steplength,int((muend-muuse)/steplength)+1)
            if muuse!=1000:
                startingdict={
                    "range":[mu1,muuse],
                    "T":float(T),
                    "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                }
                self.starting_point_list.append(startingdict)
    
    def compute_phase_diagram(self,Tstart,mustart,muend,steplength=0.1):
        #first scan phase boundary starting point
        self.scan_phase_boundary(Tstart,mustart,muend,steplength)
        print("scan finished, now start search")  
        while len(self.starting_point_list)>0:    #dictionary list of starting points
            (new_phase_boundary,signal)=self.trace_phase_boundary(self.starting_point_list[0])
            print("search finished, will move on to next starting point possible")
            self.phase_boundary_list.append(new_phase_boundary)
            self.starting_point_list.pop(0)       #after each run remove the starting point used
    
    def plot_phase_diagram0(self,filename):
        for i in range(len(self.phase_boundary_list)):
            listuse=self.phase_boundary_list[i]
            plt.plot(listuse[0],listuse[2])
            plt.plot(listuse[1],listuse[2])
        plt.savefig(filename)

    def output_optimization(self,result,T,component_comp):
        brute=0
        if brute:
            svm=self.compute_site_variable_matrix(result,component_comp,T) 
            print(" site variable is "+str(svm[1])+" composition is "+str(component_comp))
        else:
            svm=self.compute_site_variable_matrix(result.x,component_comp,T) 
            print("free energy is "+str(result.fun)+" site variable is "+str(svm[1])+" composition is "+str(component_comp))

    def output_potential(self,result,T,component_comp):
        svm=np.zeros(self.clustersize)
        for i in range(self.clustersize):
            svm[i]=np.exp(result.x[i]/T)
        #print("free energy is "+str(result.fun)+" site variable is "+str(svm[1])+" at temperature "+str(T))
        print("free energy is "+str(result.fun+mu*component_comp[1])+" site variable is "+str(svm)+" composition is "+str(component_comp)+" potential is "+str(mu))

    def print_output(self):
        print("hello world")

if __name__ == '__main__':
    print("FYLCVM code")
    E=[0,-3,-4,-3,0]
    CVM=FYLCVM(2,4,E,1)
    print(CVM.shape)
    length=20
    length+=1
    space="mu"
    space='skip'
    if space=="mu":
        #T=1.0
        #result=CVM.optimize_grand_potential(1.0,7.6,"BFGS")        
        #(total_energy,composition,F,E)=CVM.compute_grand_potential_output(result.x,T,7.5)
        #CVM.output_optimization(result,T,composition)
        CVM.trace_phase_boundary(1.2,1.5,0.02,7.0)
        #(mustart,muend,x1,x2,E1,E2,result1,result2)=CVM.search_phase_boundary(1.0,7.0,0.1,10)
        #print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2))
        
        '''muspace=np.linspace(7.5,7.6,length,endpoint=True)
        xspace=np.zeros(length)
        count=0
        for mu in muspace:
            result=CVM.optimize_grand_potential(T,mu,"NM")
            (potential,composition,E)=CVM.compute_grand_potential_output(result.x,T,mu)
            print("energy is "+str(E))
            xspace[count]=composition[0]
            CVM.output_potential(result,T,composition)
            #if count>1:
                #print((xspace[count]-xspace[count-1])/(xspace[count-1]-xspace[count-2]))
            count+=1
        plt.plot(muspace,xspace,marker="o")
        plt.show()'''
    
    elif space=="x":
        Fspace=np.zeros(length)
        Sspace=np.zeros(length)
        muspace=np.zeros(length)
        spspace=np.zeros(length)
        Cspace=np.linspace(0.36,0.41,length,endpoint=True)
        count=0
        T=1
        for i in Cspace:
            composition=np.array([i,1-i])
            result=CVM.optimize_free_energy(T,composition,"NM")
            CVM.output_optimization(result,T,composition)
            Fspace[count]=result.fun
            Cspace[count]=i
            (F,S,mu)=CVM.compute_total_energy_output(result.x,T,composition)
            Sspace[count]=S
            muspace[count]=mu
            spspace[count]=result.x[0]
            count+=1
            #print("entropy is "+str(S))

        plt.plot(Cspace,Fspace,marker="o")
        #plt.show()

    elif space=="T":
        composition=np.array([0.15,0.85])
        count=0
        Tspace=np.linspace(1.0,2.0,41)
        Fspace=np.zeros(41)
        Sspace=np.zeros(41)
        for T in Tspace:
            result=CVM.optimize_free_energy(T,composition)
            CVM.output_optimization(result,T,composition)
            (F,S)=CVM.compute_total_energy_output(result.x,T,composition)
            Sspace[count]=S  
            count+=1  
            #print("free energy is "+str(F)+"entropy is "+str(S)+" at T="+str(T))

        plt.plot(Tspace,Sspace,marker="o")
        plt.show()
    
    elif space=="point":
        composition=np.array([0.38,0.62])
        T=1
        result=CVM.optimize_free_energy(T,composition,"brute")
        #print(result)
        CVM.output_optimization(result,T,composition)
        #(F,S,mu)=CVM.compute_total_energy_output(result.x,T,composition)


    '''
    #result=CVM.optimize_free_energy(1.0,composition)
    composition=np.array([0.4,0.6])
    T=1
    F=1000
    best_site_variable=np.zeros(3)
    for i in range(400):
        for j in range(300):
            site_variable_input=np.array([1+i*0.01,1+i*0.01,0.01+j*0.01])
            Fnew=CVM.compute_total_energy(site_variable_input,T,composition)
            if Fnew<F:
                F=Fnew
                best_site_variable=site_variable_input
    svm=CVM.compute_site_variable_matrix(best_site_variable,composition,T)
    print("F is "+str(F))
    print(svm)
    
    Tspace=np.linspace(1.0,2.5,61)
    for T in Tspace:
        svm=CVM.compute_site_variable_matrix(np.array([np.log(6.95524423),np.log(6.95524423),np.log(6.95524423)]),composition,T)
        print("T is "+str(T)+"  "+str(svm[1][3]))
    
    #test first with fixed mu
    '''
    #free energy is -7.251473550123585 site variable is [8.04117739 8.04117739 0.17299931 8.04118063] composition is [0.38 0.62]