#!/usr/bin/python

import numpy as np
import os
import sys
import itertools
import utility
import plotter
import mapping
import matplotlib.pyplot as plt
import freeenergy           #This is cython library
from scipy.optimize import minimize
from scipy.optimize import brute
from scipy.optimize import basinhopping
from multiprocessing import Pool
from matplotlib import cm
from fylcvm import FYLCVM
#svm stands for site_variable_matrix

class FCC(FYLCVM):       #sub class for FCC
    def __init__(self,inputs,control_dict):                #initialize all variables
        super().__init__(inputs,control_dict)
        self.lattice="FCC"
        self.map_basic_cluster_energy(inputs.basic_cluster_energy)
        self.local_para=inputs.local_energy_parameter

    def map_basic_cluster_energy(self,E):     #map the basic cluster energy to high dimensional array
        self.basicclusterenergy=np.zeros(self.shape)                  #use high dimensional array, input by hand now
        if self.clustersize==4:
            for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
                self.basicclusterenergy[i][j][k][l]=E[i+k+j+l]
    
    def map_vibrational_energy(self,T):        #Assume symmetry wAA=wBB
        Fmat=[0,1.5*T*self.R*np.log(self.vib_para*self.local_para[1]),2*T*self.R*np.log(self.vib_para*self.local_para[2])
              ,1.5*T*self.R*np.log(self.vib_para*self.local_para[3]),0]
        #Fmat=[0,1.5*T*self.R*np.log(self.vib_para),2*T*self.R*np.log(self.vib_para*1.05),1.5*T*self.R*np.log(self.vib_para),0]
        for i,j,k,l in itertools.product(range(self.component), repeat=4):       #eith use an elegant way or write everything without loop at all
            self.vibrationalenergy[i][j][k][l]=Fmat[i+k+j+l]
    
    def map_real_energy(self,T):
        print("Map energy from energy files")     #just a place holder
        Fclass=mapping.F_FCC("FCC","energy.in")
        Emat=[0,Fclass.AAAB(T),Fclass.AABB(T),Fclass.ABBB(T),0]
        for i,j,k,l in itertools.product(range(self.component), repeat=4): 
            self.totalenergy[i][j][k][l]=Emat[i+k+j+l]

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
    
    def compute_partition_function_real(self,spm,T):
        partition=0
        for i,j,k,l in itertools.product(range(self.component), repeat=4):
            partition+=np.exp((spm[i][0]+spm[j][1]+spm[k][2]+spm[l][3]-self.totalenergy[i][j][k][l])/(self.R*T))
        return partition

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

        two_body_energy=freeenergy.compute_binary_cluster_energy(prob,"FCC")*self.R*T

        #point cluster  
        pointenergy,point_prob=freeenergy.compute_point_cluster_energy(prob,"FCC")
        pointenergy=pointenergy*self.R*T
        #point_prob=np.zeros((self.component,self.clustersize))

        composition=np.zeros(self.component)
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                composition[j]+=point_prob[j][i]*0.25
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=self.R*T*np.log(partition)
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        
        total_energy=2*basic_cluster_energy+two_body_energy+pointenergy-composition[1]*mu+elastic_energy
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

        two_body_energy=freeenergy.compute_binary_cluster_energy(prob,"FCC")*self.R*T

        #point cluster  
        pointenergy,point_prob=freeenergy.compute_point_cluster_energy(prob,"FCC")
        pointenergy=pointenergy*self.R*T
        #basic cluster start with tetra as default    
        
        composition=np.zeros(self.component)    
        basic_cluster_energy=0                           
        for i in range(self.clustersize):
            for j in range(self.component):
                basic_cluster_energy+=point_prob[j][i]*spm[j][i]
        basic_cluster_energy-=T*self.R*np.log(partition)

        E=0
        #for i,j,k,l in itertools.product(range(self.component), repeat=4):
        #    E+=self.basicclusterenergy[i][j][k][l]*prob[i][j][k][l]
        
        elastic_energy=self.elastic_parameter*composition[0]*composition[1]
        grand_potential=2*basic_cluster_energy-two_body_energy+1.25*pointenergy-composition[1]*mu+elastic_energy
        F=grand_potential+composition[1]*mu
        #entropy=(2*np.sum(self.basicclusterenergy*prob)-F)/T
        return grand_potential,composition,F,E
    
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
            T=round(self.dT*count*phb.direction+phb.Tstart,2)
            print("T is "+str(T))
            if T<self.Tmin:
                break
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
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-phb.direction)
                        phb.status="exit"
                        phb.signal=1
                        return phb        
                if self.identify_phase(result2.x)!=highphase:                  #new phase found
                    print("new phase "+self.identify_phase(result2.x)+" detected, follow two new boundary")
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,muend,self.dmurough,10)
                    if muend:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                        phb.status="exit"
                        phb.signal=2
                    else:
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,muendsave,self.dmurough,20)
                        if mustart:
                            self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-phb.direction)
                            phb.status="exit"
                            phb.signal=2
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                            phb.signal=-1
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
                    print("can this even happen?")
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,20)
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                        phb.signal=-1
                        return phb
                    else:
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-phb.direction)
                        phb.status="exit"
                        phb.signal=1
                        return phb   
                if self.identify_phase(result1.x)!=lowphase:
                    print("new phase "+self.identify_phase(result1.x)+" detected, follow two new boundary")
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustart,-self.dmurough,10)
                    if muend:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                        phb.status="exit"
                        phb.signal=2
                    else:
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,10)
                        if mustart:
                            self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-phb.direction)
                            phb.status="exit"
                            phb.signal=2
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                            phb.signal=-1
                    return phb
            else:   
                mustartuse=mustart       
                print("new phase "+currentphase+" detected, start two new search")
                (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,self.dmurough,20)    #one forward and one reverse
                #print("mustart is "+str(mustart))
                #print("muend is "+str(muend))
                if muend and muend!=1000:   #edit later for unlikely boundary condition
                    self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,-self.dmurough,20)
                    if muend and muend!=1000:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                        phb.status="exit"
                        phb.signal=2
                    elif muend==1000:
                        phb.status="exit"
                        phb.signal=1
                    elif muend==False:
                        print("we are here?")
                        (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,mustartsave,-self.dmurough,10)
                        if mustart:
                            self.start_new_phb(self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,phb.direction)
                            phb.status="exit"
                            phb.signal=2
                        else:
                            print("phase boundary end unexpectedly")
                            phb.status="error"
                            phb.signal=-1
                elif muend==False:  #we are here now
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T,mustartuse,-self.dmurough,20)    #The reverse search must have boundary
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                        phb.signal=-1
                    else:
                        self.start_new_phb(T,mustart,muend,result1,result2,x1,x2,phb.direction)
                        phb.status="exit"
                    #This line has error too
                    (mustart,muend,x1,x2,E1,E2,result1,result2)=self.search_phase_boundary_v1(T-self.dT*phb.direction,muendsave,self.dmurough,20)    #The reverse search must have boundary
                    if not mustart:
                        print("phase boundary end unexpectedly")
                        phb.status="error"
                        phb.signal=-1
                    else:
                        self.start_new_phb(T-self.dT*phb.direction,mustart,muend,result1,result2,x1,x2,-phb.direction)
                        phb.status="exit"
                        phb.signal=2
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
        reverse=0
        if mustart>muend:
            reverse=1
            mustart,muend=muend,mustart         #make sure mustart is always smaller
            phasestart,phaseend=phaseend,phasestart
        resultlist=[]
        while np.abs(mustart-muend)>self.mutol:
            mumid=0.5*(mustart+muend)
            result=self.optimize_grand_potential(T,mumid,"BFGS_v2") 
            resultlist.append(result)
            phasemid=self.identify_phase(result.x)
            if phasemid!=phasestart and phasemid!=phaseend:        #need recursion
                print("new phase "+phasemid+" detected? Steplength for phb search might be too big")
                if reverse:
                    return self.search_phase_boundary_bisection(T,mumid,muend,phasemid,phaseend)
                else:
                    return self.search_phase_boundary_bisection(T,mustart,mumid,phasestart,phasemid)
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
                    self.phb_status_list.append(self.node_list[i].signal)
        f=open(phbfile_name,"w")
        print(self.phase_boundary_list,file=f)
        f.close()
        utility.replace_word_in_file(phbfile_name,"array","np.array")          #make the output file that can be directly used by plotter
        #print(self.phase_boundary_list)
        #signal_list=[2,0,0,0]
        print(self.phb_status_list)
        #plotter.plot_phase_diagram(self.phase_boundary_list,self.phb_status_list,self.dT,Tstart,phase_diagram_name)

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

    def plot_point_region(self,Trange,murange,Tstep,mustep):
        data = np.zeros((0, 3))
        muTdata = np.zeros((0, 3))
        Tarray=np.arange(Trange[0],Trange[1],Tstep)
        muarray=np.arange(murange[0],murange[1],mustep)
        for T in Tarray:
            for mu in muarray:
                result=self.optimize_grand_potential(T,mu,"BFGS_v2")
                GP,composition,F,E=self.compute_grand_potential_output(result.x,T,mu)
                phase=self.identify_phase(result.x)
                data=np.append(data,np.array([[composition[0], T, phase]]),axis=0)
                muTdata=np.append(muTdata,np.array([[mu, T, phase]]),axis=0)
        #print(data)
        plotter.plot_area(data)
        plotter.plot_area(muTdata,True)

    def plotGx(self,T,mustart,muend,dmu=0.05):
        muarr=np.arange(mustart,muend,dmu)
        Farray=np.array([])
        xarray=np.array([])
        for mu in muarr:
            result=self.optimize_grand_potential(T,mu,"BFGS_v2")
            GP,composition,F,E=self.compute_grand_potential_output(result.x,T,mu)
            Farray=np.append(Farray,F)
            xarray=np.append(xarray,composition[0])
        plt.scatter(xarray,Farray)
        plt.show()
        f=open("gx.txt","w")
        print(np.vstack((xarray,Farray)),file=f)
        f.close()
if __name__ == '__main__':
    print("you are in FYLCVM code?")