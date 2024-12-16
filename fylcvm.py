#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import plotter
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import brute
from scipy.optimize import basinhopping
from multiprocessing import Pool
from matplotlib import cm
#svm stands for site_variable_matrix

class FYLCVM:       #base class for low symmetry structure
    def __init__(self,component,maxsize,E,vibration_parameter,local_parameter,elastic_para,control_dict):                #initialize all variables
        #print("Initialize FYL-CVM\n")
        kB=1.380649e-23  #boltzman constant, unit:J/K
        N=6.02e23;       #number of particles, the size of the system, here we use 1 mole
        self.lattice="monoclinic"
        self.component=component
        self.clustersize=maxsize
        self.shape=[self.component]*self.clustersize    #ignore using basic cluster with different size first
        self.map_basic_cluster_energy(E)
        self.vibrationalenergy=np.zeros(self.shape)
        self.totalenergy=np.zeros(self.shape)
        self.vib_para=vibration_parameter             #vibration parameter with default of 1
        self.local_para=local_parameter
        self.elastic_parameter=elastic_para
        self.starting_point_list=[]
        self.node_list=[]
        self.phase_boundary_list=[]
        self.mu_Tlist=[]
        self.phb_status_list=[]
        self.R=control_dict['R']
        self.dmurough=control_dict['dmurough']
        self.dmuprecise=control_dict['dmuprecise']
        self.dmuscan=control_dict['dmuscan']
        self.dT=control_dict['dT']
        self.Tmin=control_dict['Tstart']
        self.Tmax=control_dict['Tmax']
        self.mutol=control_dict['mutol']
        self.mumin=control_dict['mumin']
        self.mumax=control_dict['mumax']
    
    def update_parameter(self,keyword,value,position=1):
        if keyword=="e":          #elastic
            self.elastic_parameter=value
        elif keyword=="vib":
            self.vib_para=value
        elif keyword=="local":
            self.local_para[position]=value
        else:
            print("no parameter to update")

    def map_basic_cluster_energy(self,E):     #map the basic cluster energy to high dimensional array
        self.basicclusterenergy=np.zeros(self.shape)                  #use high dimensional array, input by hand now
        if self.clustersize==4:
            for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
                self.basicclusterenergy[i][j][k][l]=E[i+k+j+l]
    
    '''def map_vibrational_energy(self,T):        #w is w matrix normalized
        wmat=self.vib_matrix
        pab=1.05
        Fmat=[-3*T*np.log(T/wmat[0]),-1.5*T*np.log(T/wmat[0])-1.5*T*np.log(T/wmat[1]),-0.5*T*np.log(T/wmat[0])
             -2*T*np.log(T/wmat[1])-0.5*T*np.log(T/wmat[2]),-1.5*T*np.log(T/wmat[1])-1.5*T*np.log(T/wmat[2]),-3*T*np.log(T/wmat[2])]
        for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
            self.vibrationalenergy[i][j][k][l]=Fmat[i+k+j+l]'''
    
    def map_real_energy(self,T):
        print("This depends on crystal structure")     #just a place holder
    
    def map_vibrational_energy(self,T):        #Assume symmetry wAA=wBB
        Fmat=[0,1.5*T*self.R*np.log(self.vib_para*self.local_para[1]),2*T*self.R*np.log(self.vib_para*self.local_para[2])
              ,1.5*T*self.R*np.log(self.vib_para*self.local_para[3]),0]
        #Fmat=[0,1.5*T*self.R*np.log(self.vib_para),2*T*self.R*np.log(self.vib_para*1.05),1.5*T*self.R*np.log(self.vib_para),0]
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
            nominator+=(composition[1]-(1/self.clustersize)*(i+j+k+0))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][0]/(self.R*T))
            denominator+=(-composition[1]+(1/self.clustersize)*(i+j+k+1))*svm[i][0]*svm[j][1]*svm[k][2]*np.exp(-self.basicclusterenergy[i][j][k][1]/(self.R*T))
        svm[1][3]=nominator/denominator
        return svm
        
        #merge vibrational energy into basic cluster energy
    def compute_partition_function(self,svm,T,Fvib):
        if svm[0][0]==1:                 #matrix is for site variable
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp((-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))
            return partition
        else:                            #matrix is for site potential, svm is actually spm
            partition=0
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                partition+=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))
            return partition
        
    def compute_basic_cluster_prob(self,partition,svm,T,Fvib):        #svm is site potential matrix
        basic_cluster_prob=np.zeros(self.shape)
        if svm[0][0]==1:                 #matrix is for site variable
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                basic_cluster_prob[i][j][k][l]=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]*np.exp((-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))/partition
        else:                            #matrix is for site potential, svm is actually spm
            for i,j,k,l in itertools.product(range(self.component), repeat=4):
                basic_cluster_prob[i][j][k][l]=np.exp((svm[i][0]+svm[j][1]+svm[k][2]+svm[l][3]-self.basicclusterenergy[i][j][k][l]-Fvib[i][j][k][l])/(self.R*T))/partition
        return basic_cluster_prob
    
    def compute_additional_energy(self,vibrational,elastic,electronic):     #vibrational contains the ratio
        if vibrational:
            print("helloworld")

    def identify_phase(self,site_potential_input):                           #identify which phase is showing up
        count=0
        for i in range(self.clustersize):
            if site_potential_input[i]==0:
                return "A1"
            for j in range(i+1,self.clustersize):
                if np.abs((site_potential_input[i]-site_potential_input[j])/site_potential_input[i])<3.0e-2 or np.abs(site_potential_input[i]-site_potential_input[j])<0.01:
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
            #print("negative partition function")
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
                    two_body_energy+=probij*np.log(probij)*self.R*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*self.R*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=self.R*T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=self.R*T*np.log(partition)
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        return total_energy

    def compute_total_energy_brute(self,site_variable_input,T,component_comp):    #reduce computational cost for brute method using special function
        site_variable_input=np.append(site_variable_input,site_variable_input[0])
        svm=self.compute_site_variable_matrix(site_variable_input,component_comp,T) 
        partition=self.compute_partition_function(svm,T,self.vibrationalenergy)
        prob=self.compute_basic_cluster_prob(partition,svm,T,self.vibrationalenergy)

        if partition<0:
            #print("negative partition function")
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
                    two_body_energy+=probij*np.log(probij)*self.R*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*self.R*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=self.R*T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=self.R*T*np.log(partition)
        
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
                    two_body_energy+=probij*np.log(probij)*self.R*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*self.R*T
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=self.R*T*(point_prob[j][i]*np.log(svm[j][i]))
        basic_cluster_energy-=self.R*T*np.log(partition)
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
        mu=0
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            product=svm[i][0]*svm[j][1]*svm[k][2]*svm[l][3]
            mu+=product*np.log(product)*np.exp(self.basicclusterenergy[i][j][k][l])*T

        return total_energy,entropy,mu/partition

    def compute_grand_potential(self,site_potential_input,T,mu):    #compute the grand cononical potential
        self.map_vibrational_energy(T)
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
                    two_body_energy+=probij*np.log(probij)*self.R*T
        #point cluster  
        pointenergy=0
        point_prob=np.zeros((self.component,self.clustersize))
        composition=np.zeros(self.component)
        for i in range(self.clustersize):
            for j in range(self.component):
                position_type_matrix=np.zeros((2,1))
                position_type_matrix[0][0]=i
                position_type_matrix[1][0]=j
                point_prob[j][i]=utility.get_cluster_prob(prob,position_type_matrix)
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*self.R*T
                composition[j]+=point_prob[j][i]*0.25
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=self.R*T*np.log(partition)
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy-composition[1]*mu+elastic_energy
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
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
                    two_body_energy+=probij*np.log(probij)*T*self.R
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
                pointenergy+=point_prob[j][i]*np.log(point_prob[j][i])*T*self.R
                composition[j]+=point_prob[j][i]*0.25
        #basic cluster start with tetra as default        
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=T*self.R*np.log(partition)

        E=0
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            E+=self.basicclusterenergy[i][j][k][l]*prob[i][j][k][l]
        
        total_energy=2*basic_cluster_energy-two_body_energy+1.25*pointenergy-composition[1]*mu
        F=total_energy-composition[1]*mu
        entropy=(2*np.sum(self.basicclusterenergy*prob)-total_energy)/T
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
        elif method=="BFGS":
            A1_1=np.array([0,0,0,0])
            A1_2=np.array([20,20,20,20])
            A1_3=np.array([-20,-20,-20,-20])
            L12_1=np.array([-10,-10,-10,10])
            L12_2=np.array([3,3,3,-1.6])
            L12_3=np.array([10,10,10,-10])
            L10_1=np.array([15,15,-15,-15])   
            L10_2=np.array([-15,-15,15,15])  
            L10_2=np.array([5,5,-5,-5])
            L10_2=np.array([-5,-5,5,5]) 
            guesslist=(A1_1,A1_2,A1_3,L12_1,L12_2,L12_3,L10_1,L10_2)
            #result=minimize(self.compute_grand_potential,guess,method='BFGS',args=(T,mu))
            result=self.multiBFGS(guesslist,(T,mu))
        elif method=="BFGS_v2":
            '''guess_1=np.array([5,-5])
            guess_2=np.array([-5,5])
            L10_1=np.array([15,-15,])   
            L10_2=np.array([-15,15,])'''
            guess_1=np.array([2,-1])
            guess_2=np.array([-1,2])
            L10_1=np.array([5,-5,])   
            L10_2=np.array([-5,5,])
            guesslist=(guess_1,guess_2,L10_1,L10_2)  
            result,phase=self.multiBFGS_v2(guesslist,(T,mu))
            result.x=utility.map_variable_to_phase(result.x,phase,self.clustersize)         #change later if everything works perfectly, return result.x with size 4
            #print("at potential mu="+str(mu)+" best phase is "+self.identify_phase(result.x))
        return result
    
    def plot_potential_surface(self,T,mu,phase):
        print("plot the potential surface based on brute force plot result")
        if phase=="A1":
            boundary=[slice(-10,10,0.02)]
        else:
            boundary=(slice(-50,50,0.2),slice(-50,50,0.2))
        result=brute(self.brute_minimizer,boundary,args=(T,mu,phase),workers=-1,full_output=True,finish=None)
        print("minimal is at "+str(result[0])+" minimal value is "+str(result[1]))
        if phase=="A1":
            fig, ax = plt.subplots(subplot_kw={})
            ax.plot(result[2],result[3])
        else:
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            surf = ax.plot_surface((result[2])[0],(result[2])[1],result[3],linewidth=0, cmap=cm.coolwarm,antialiased=False)
        ax.set_title("potential surface with mu="+str(mu)+" at T="+str(T)+" for phase "+phase)
        fig.savefig("mu="+str(mu)+"_"+phase)
    
    def brute_minimizer(self,input_variable,T,mu,phase):
        site_potential_input=np.zeros(self.clustersize)
        if phase=="A1":
            for i in range(self.clustersize):site_potential_input[i]=input_variable[0] 
        elif phase=="L12":
            site_potential_input[0]=site_potential_input[1]=site_potential_input[2]=input_variable[0]
            site_potential_input[3]=input_variable[1]
        elif phase=="L10":
            site_potential_input[0]=site_potential_input[1]=input_variable[0]
            site_potential_input[2]=site_potential_input[3]=input_variable[1]
        return self.compute_grand_potential(site_potential_input,T,mu)
    
    def multiBFGS(self,guesslist,args,space="mu"):             #input two tuple
        Fmin=1000
        argslist=[args]*len(guesslist)
        inputlist=zip(guesslist,argslist)
        if space=="mu":
            with Pool(processes=8) as pool:
                resultlist=pool.starmap(self.minimizermu,inputlist)
                pool.close()
        else:
            with Pool(processes=8) as pool:
                resultlist=pool.starmap(self.minimizerx,inputlist)
                pool.close()

        for result in resultlist:
            #print("free energy is "+str(result.fun)+" site variable is "+str(result.x))
            if result.fun<Fmin:
                if space!="mu" or self.identify_phase(result.x)!="L12+L10":                   #ignore L12+L10 for now
                    Fmin=result.fun
                    bestresult=result
        return bestresult
    
    def multiBFGS_v2(self,guesslist,args):             #input two tuple
        Fmin=10
        argsL12=args+("L12",)
        argsL10=args+("L10",)
        argslist=(argsL12,argsL12,argsL12,argsL12,argsL10,argsL10,argsL10,argsL10)
        #print(argslist)
        #longguesslist=guesslist*2
        #print(longguesslist)
        inputlist=zip(guesslist*2,argslist)

        with Pool(processes=8) as pool:
            resultlist=pool.starmap(self.minimizerv2,inputlist)
            pool.close()

        for result in resultlist:
            #print("phase is "+result[1])
            if result[0].fun<Fmin:
                Fmin=result[0].fun
                bestresult=result[0]
                phase=result[1]
        
        return bestresult,phase
    
    def minimizermu(self,initialvalue,args,method='BFGS'):      #use a minimizer that can be called by pool.starmap
        result=minimize(self.compute_grand_potential,initialvalue,method=method,args=args)
        return result
    
    def minimizerv2_L12(self,initialvalue,args,method='BFGS'):      #reduce number of variables
        result=minimize(self.potential_computerv2_L12,initialvalue,method=method,args=args)
        return result

    def minimizerv2(self,initialvalue,args,method='BFGS'):      #return phase
        result=minimize(self.potential_computerv2,initialvalue,method=method,args=args)
        return result,args[2]
    
    def minimizerv2_L10(self,initialvalue,args,method='BFGS'):      #reduce number of variables
        print("check if minimizer is runned")
        result=minimize(self.potential_computerv2_L10,initialvalue,method=method,args=args)
        return result

    def potential_computerv2(self,initialvalue,T,mu,phase):
        sitepotential_input=np.zeros(self.clustersize)
        if phase=="L12":
            sitepotential_input[0]=sitepotential_input[1]=sitepotential_input[2]=initialvalue[0]
            sitepotential_input[3]=initialvalue[1]
        elif phase=="L10":
            sitepotential_input[0]=sitepotential_input[1]=initialvalue[0]
            sitepotential_input[2]=sitepotential_input[3]=initialvalue[1]
        return self.compute_grand_potential(sitepotential_input,T,mu)
        
    def potential_computerv2_L12(self,initialvalue,T,mu):
        sitepotential_input=np.zeros(self.clustersize)
        sitepotential_input[0]=sitepotential_input[1]=sitepotential_input[2]=initialvalue[0]
        sitepotential_input[3]=initialvalue[1]
        return self.compute_grand_potential(sitepotential_input,T,mu)

    def potential_computerv2_L10(self,initialvalue,T,mu):    
        sitepotential_input=np.zeros(self.clustersize)
        sitepotential_input[0]=sitepotential_input[1]=initialvalue[0]
        sitepotential_input[2]=sitepotential_input[3]=initialvalue[1]
        return self.compute_grand_potential(sitepotential_input,T,mu)

    def minimizerx(self,initialvalue,args,method='BFGS'):      #use a minimizer that can be called by pool.starmap
        result=minimize(self.compute_total_energy,initialvalue,method=method,args=args)
        return result   

    def optimize_free_energy(self,T,component_comp,method):         #main optimization function, optimize site variables at given T and composition
        guess=np.array([2.26,2.26,-1.84])
        #initial_guess=np.ones((self.clustersize-1)*(self.component-1))
        if method=="NM":
            initial_guess2=np.array([[5,5,-2],[5,5,-1],[4,4,0],[4,4,-1]])
            options={
                'initial_simplex':initial_guess2
            }
            positive=((-10,10),(-10,10),(-10,10))
            #result=minimize(CVM.compute_total_energy,guess,method='Nelder-Mead',args=(T,component_comp),bounds=positive,tol=1.0e-6)
            result=minimize(self.compute_total_energy,guess,method='Nelder-Mead',options=options,args=(T,component_comp),bounds=positive,tol=1.0e-6)
        elif method=="basinhopping":
            minimizer_kwargs={"args":(T,component_comp)}
            result=basinhopping(self.compute_total_energy,guess,minimizer_kwargs=minimizer_kwargs)
        elif method=="brute":
            boundary=(slice(-10,10,0.02),slice(-10,10,0.02))
            result=brute(self.compute_total_energy_brute,boundary,args=(T,component_comp),workers=-1,full_output=False)
        elif method=="BFGS":
            A1_1=np.array([2.5,2.5,2.5])
            A1_2=np.array([0,0,0])
            L12_1=np.array([3.2,3.2,-1])
            L10_1=np.array([3.5,3.5,-0.5])   
            guesslist=(A1_1,A1_2,L12_1,L10_1)
            #result=minimize(self.compute_grand_potential,guess,method='BFGS',args=(T,mu))
            result=self.multiBFGS(guesslist,(T,component_comp),"x")

        return result
    
    def find_phase_boundary(self,Tstart,Tstop,composition,method="BFGS"):
        T=Tstart
        count=0
        while T<Tstop:
            result=self.optimize_free_energy(T,composition,method)
            self.output_optimization(result,T,composition,method)
            T+=self.dT
    
    #trace_phase_boundary returns the phase boundary between two specific phases
    def trace_phase_boundary(self,startingdict):           #always move in the direction of incresing mu
        lowphase=startingdict['phase'][0]
        highphase=startingdict['phase'][1]
        Tmin=startingdict['T']
        mustart=startingdict['range'][0]
        #if mustart==False:
        #    print("invalid input")
        #    return False,False,""
        print("trace phase boundary between "+lowphase+" and "+highphase)
        x1mat=np.array([])
        x2mat=np.array([])
        Tspace=np.array([])               #flexable value
        muspace=np.array([])  
        T=Tmin
        count=0

        while (T<=self.Tmax):
            T=self.dT*count+Tmin
            mustart+=self.compute_dmu()             #edit after compute_dmu() has value
            result=self.optimize_grand_potential(T,mustart,"BFGS_v2")
            currentphase=self.identify_phase(result.x)
            #print("current phase is "+currentphase)
            if currentphase==lowphase:
                #print("lowphase")
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,self.dmuprecise,150)
                if muend==False:
                    print("phase boundary ends")
                    return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),0
                if self.identify_phase(result2.x)!=highphase:
                    print("new phase detected, follow two new boundary")
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(newstartingdict)
                    print("new starting dict is")
                    print(newstartingdict)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,muend,self.dmurough,100)
                    if muend:
                        newstartingdict={
                            "range":[mustart,muend],
                            "T":T,
                            "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                        }
                        print("new starting dict is")
                        print(newstartingdict)
                        self.starting_point_list.append(newstartingdict)
                    else:
                        print("phase boundary might go downward")
                    return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),2
            elif currentphase==highphase:
                #print("highphase")
                #search in reverse so that result is also in reverse
                (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-self.dmuprecise,150)
                if muend==False:
                    print("phase boundary ends")
                    return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),0
                if self.identify_phase(result1.x)!=lowphase:
                    print("new phase detected, follow two new boundary")
                    newstartingdict={
                        "range":[mustart,muend],
                        "T":T,
                        "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    print("new starting dict is")
                    print(newstartingdict)
                    self.starting_point_list.append(newstartingdict)
                    (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-self.dmurough,100)
                    if muend:
                        newstartingdict={
                            "range":[mustart,muend],
                            "T":T,
                            "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                        }
                        print("new starting dict is")
                        print(newstartingdict)
                        self.starting_point_list.append(newstartingdict)
                    else:
                        print("new phase boundary goes downward")
                    return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),2
            else:
                print("new phase detected, start two new search")
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,mustart,self.dmurough,100)    #one forward and one reverse
                if muend:
                    newstartingdict={
                            "range":[mustart,muend],
                            "T":T,
                            "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                        }
                    print("new starting dict is")
                    print(newstartingdict)
                    self.starting_point_list.append(newstartingdict)
                (muend,mustart,x2,x1,E2,E1,result2,result1)=self.search_phase_boundary(T,mustart,-self.dmurough,100)
                if muend:
                    newstartingdict={
                            "range":[mustart,muend],
                            "T":T,
                            "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                        }
                    print("new starting dict is")
                    print(newstartingdict)
                    self.starting_point_list.append(newstartingdict)
                return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),2
            #print("mustart is "+str(mustart)+" muend is "+str(muend)+" x1 is "+str(x1)+" x2 is "+str(x2)+" at T="+str(T))
            print("mu is between {:.2f} and {:.2f} composition is between {:.3f} and {:.3f} at T={:.2f}".format(mustart,muend,x1,x2,T))
            #print("sitepot1 is "+str(result1.x)+" sitepot2 is "+str(result2.x))     
            x1mat=np.append(x1mat,x1)
            x2mat=np.append(x2mat,x2)
            Tspace=np.append(Tspace,T)
            muspace=np.append(muspace,0.5*(mustart+muend))
            count+=1
        
        return np.stack((x1mat,x2mat,Tspace)),np.stack((muspace,Tspace)),1

    def trace_phase_boundary_v1(self,startingdict):
        print("helloworld placeholder")

    def compute_dmu(self):
        return 0

    def search_phase_boundary(self,T,mustart,steplength,steps):
        if steps>1000:            #max of 1000 steps
            steps=1000
        cdata=np.zeros(steps)
        for i in range(steps):
            muuse=mustart+i*steplength
            result=self.optimize_grand_potential(T,muuse,"BFGS_v2")  
            #print("current phase is "+str(self.identify_phase(result.x))+" at potential "+str(muuse)+" para is "+str(result.x))
            (potential,composition,F,E)=self.compute_grand_potential_output(result.x,T,muuse)
            cdata[i]=composition[0]
            currentphase=self.identify_phase(result.x)
            if i>0:
                if currentphase!=phasesave:
                    #print("phase transition between "+currentphase+" and "+phasesave+"\nmu from "+str(muuse-steplength)+" to "+str(muuse))
                    print("phase transition between {} and {}".format(phasesave,currentphase)) 
                    #print("site potential is "+str(result.x)+" and "+str(resultsave.x))             
                    return muuse-steplength,muuse,cdata[i-1],cdata[i],Esave,E,resultsave,result
            phasesave=currentphase
            Esave=E
            resultsave=result
        print("can't find phase boundary")
        return False,False,cdata[i-1],cdata[i],Esave,E,resultsave,result
    
    def scan_phase_boundary(self,T,mustart,muend,scan=True):        #scan the phase boundary across a wide 
        print("scan the phase boundary between mu={:.2f} and {:.2f} at T={:.2f}".format(mustart,muend,T))
        phblist=np.empty((0,2))
        mulist=np.empty(0)
        muuse=mustart
        while muuse and muuse<muend:                             #continue search with new muuse until boundary is reached
            (mu1,muuse,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary(T,muuse,self.dmuscan,int((muend-muuse)/self.dmuscan)+1)
            if muuse and mu1:
                if scan:
                    startingdict={
                    "range":[mu1,muuse],
                    "T":float(T),
                    "phase":[self.identify_phase(result1.x),self.identify_phase(result2.x)]
                    }
                    self.starting_point_list.append(startingdict)
                else:
                    phblist=np.append(phblist,np.array([[x1,T],[x2,T]]),axis=0)
                    mulist=np.append(mulist,np.array([mu1]))
                    mulist=np.append(mulist,np.array([muuse]))
        if scan==False and np.size(mulist)>0:                                #This is for scanning the whole phase diagram
            return phblist,np.min(mulist),np.max(mulist)
        else:
            return phblist,0,0
    
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
    
    def compute_phase_diagram(self,Tstart,mustart,muend,phase_diagram_name):
        #first scan phase boundary starting point
        placeholder=self.scan_phase_boundary(Tstart,mustart,muend)
        print("scan finished, now start search")  
        while len(self.starting_point_list)>0:    #dictionary list of starting points
            (new_phase_boundary,mu_T,signal)=self.trace_phase_boundary(self.starting_point_list[0])
            print("search finished, will move on to next starting point possible")
            self.phase_boundary_list.append(new_phase_boundary)
            self.mu_Tlist.append(mu_T)
            self.phb_status_list.append(signal)
            self.starting_point_list.pop(0)       #after each run remove the starting point used
        f=open("phase_boundary_list.txt","w")
        print(self.phase_boundary_list,file=f)
        f.close()
        utility.replace_word_in_file("phase_boundary_list.txt","array","np.array")          #make the output file that can be directly used by plotter
        plotter.plot_phase_diagram(self.phase_boundary_list,self.phb_status_list,self.dT,Tstart,phase_diagram_name)

    def plot_phase_diagram0(self,filename):
        for i in range(len(self.phase_boundary_list)):
            listuse=self.phase_boundary_list[i]
            plt.plot(listuse[0],listuse[2])
            plt.plot(listuse[1],listuse[2])
        plt.savefig(filename)

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