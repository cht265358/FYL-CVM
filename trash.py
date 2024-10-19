#!/usr/bin/python
#code trashcan

for i,j in itertools.product(range(self.component), repeat=2):
            probij=0
            for k,l in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]
            two_body_energy+=probij*np.log(probij)*T
        for i,k in itertools.product(range(self.component), repeat=2):
            probij=0
            for j,l in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]
            two_body_energy+=probij*np.log(probij)*T
        for i,l in itertools.product(range(self.component), repeat=2):
            probij=0
            for j,k in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]
            two_body_energy+=probij*np.log(probij)*T
        for j,k in itertools.product(range(self.component), repeat=2):
            probij=0
            for i,l in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]
            two_body_energy+=probij*np.log(probij)*T
        for j,l in itertools.product(range(self.component), repeat=2):
            probij=0
            for i,k in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]
            two_body_energy+=probij*np.log(probij)*T
        for k,l in itertools.product(range(self.component), repeat=2):
            probij=0
            for i,j in itertools.product(range(self.component), repeat=2):
                probij+=prob[i][j][k][l]

        #energy constraint
        penalty_parameter=1000                                #set is to large value later
        constraint_function=0
        computed_comp=np.zeros(self.component)
        for i in range(self.component):
            for j in range(self.clustersize):
                computed_comp[i]+=point_prob[i][j]/self.clustersize
            constraint_function+=penalty_parameter*np.abs((computed_comp[i])-component_comp[i])

    def search_phase_boundary(self,T,mustart,steplength,steps):
        if steps>1000:            #max of 1000 steps
            steps=1000
        cdata=np.zeros(steps)
        for i in range(steps):
            muuse=mustart+i*steplength
            result=self.optimize_grand_potential(T,muuse,"BFGS")
            (potential,composition,F,E)=self.compute_grand_potential_output(result.x,T,muuse)
            cdata[i]=composition[0]
            if i>1:
                if np.abs(((cdata[i]-cdata[i-1])/(cdata[i-1]-cdata[i-2]))-1)>0.1:
                    print("slope difference is "+str(np.abs(((cdata[i]-cdata[i-1])/(cdata[i-1]-cdata[i-2]))-1)))              
                    return muuse-steplength,muuse,cdata[i-1],cdata[i],Esave,E,resultsave,result
            Esave=E
            resultsave=result
    '''