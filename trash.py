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

'''
x2=np.array([17,8,9,8,8,8])
x3=np.array([1,2,3,4,7,6])

x=np.stack((x1,x2,x3),axis=0)
y=np.stack((x2,x3,x1),axis=0)
lists=[]
lists.append(x)
lists.append(y)
for i in range(len(lists)):
    listuse=lists[i]
    plt.plot(listuse[0],listuse[2])
    plt.plot(listuse[1],listuse[2])

plt.show()
print(x)'''

#plot brute optimize result
def compute_result(arr):
    return arr[0]**4+2.0*arr[0]*arr[1]**2+1.28*arr[1]**3

if __name__ == '__main__':
    print("test brute plot")
    boundary=(slice(-10,10,5),slice(-10,10,5))
    result=brute(compute_result,boundary,workers=-1,full_output=True)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X = np.arange(-10,10,5)
    Y = np.arange(-10,10,5)
    X, Y = np.meshgrid(X, Y)
    print((result[2])[0])
    surf = ax.plot_surface((result[2])[0],(result[2])[1],result[3],linewidth=0, antialiased=False)
    plt.show()


    def plot_phase_diagram0(self,filename):
        for i in range(len(self.phase_boundary_list)):
            listuse=self.phase_boundary_list[i]
            plt.plot(listuse[0],listuse[2])
            plt.plot(listuse[1],listuse[2])
        plt.savefig(filename)


or i in range(listlen):
        for j in range(i+1,listlen):
            if signal_list[i]!=4 and signal_list[j]!=4:
                if phase_boundary_list[i][2][-1]==phase_boundary_list[j][2][-1]:    #same end point
                    #if phase_boundary_list[i][0][0]>
                    phase_boundary_list[j][[0,1]]=phase_boundary_list[j][[1,0]]
                    list_use=np.append(phase_boundary_list[i],phase_boundary_list[j],axis=1)
                    if list_use[0][0]>list_use[0][1]:
                        list_use=np.fliplr(list_use) 
                    phb_right=UnivariateSpline(list_use[0],list_use[2],s=0.5)
                    phb_left=UnivariateSpline(list_use[1],list_use[2],s=0.5)
                    x1space=np.linspace(list_use[0][0],list_use[0][-1],1000)
                    x2space=np.linspace(list_use[1][0],list_use[1][-1],1000)
                    plt.plot(x1space,phb_right(x1space),'r')
                    plt.plot(x2space,phb_left(x2space),'r')
                    signal_list[i]=signal_list[j]=4
                    #if phase_boundary_list[i][1][-1]>phase_boundary_list[j][0][-1]:
                    #avg=0.5*(phase_boundary_list[i][1][-1]+phase_boundary_list[j][0][-1])
                    #phase_boundary_list[i]=np.append(phase_boundary_list[i],np.array([[avg],[avg],[phase_boundary_list[i][2][-1]+0.5*dT]]),axis=1)
                    #phase_boundary_list[i]=np.append(phase_boundary_list[i],phase_boundary_list[j],axis=1)
                    #phase_boundary_list[j]=np.append(phase_boundary_list[j],np.array([[avg],[avg],[phase_boundary_list[i][2][-1]]]),axis=1)
    '''
    print(phase_boundary_list)
    for i in range(listlen):
        list_use=phase_boundary_list[i]
        if len(list_use[0])>3:
            phb_right=UnivariateSpline(list_use[0],list_use[2],s=0.5)
            phb_left=UnivariateSpline(list_use[1],list_use[2],s=0.5)
            x1space=np.linspace(list_use[0][0],list_use[0][-1],1000)
            x2space=np.linspace(list_use[0][0],list_use[0][-1],1000)
            plt.plot(x1space,phb_right(x1space),'r')
            plt.plot(x2space,phb_left(x2space),'r')
        else:
            plt.plot(list_use[0],list_use[2],'r')
            plt.plot(list_use[1],list_use[2],'r')
    plt.show()
    '''

def find_intersect(mu_t1,mu_t2):
    mu_t1,flipmu1=check_positive(mu_t1)
    mu_t2,flipmu2=check_positive(mu_t2)
    fit1=UnivariateSpline(mu_t1[0],mu_t1[1])
    fit2=UnivariateSpline(mu_t2[0],mu_t2[1])
    Fmin=10
    mubest=0
    Tbest=300
    if flipmu1+flipmu2==2:    #two positive line
        ls=np.linspace(np.max(mu_t1[0][-1],mu_t2[0][-1]),np.max(mu_t1[0][-1],mu_t2[0][-1])+10,1001,endpoint=True)
    elif flipmu1==1 and flipmu2==-1:
        ls=np.linspace(mu_t1[0][-1],mu_t2[0][0],500)
    elif flipmu1==-1 and flipmu2==1:
        ls=np.linspace(mu_t2[0][-1],mu_t1[0][0],500)
    else:
        ls=np.linspace(np.min(mu_t1[0][0],mu_t2[0][0]),np.min(mu_t1[0][0],mu_t2[0][0])-10,1001,endpoint=True)
    for i in ls:
        F=np.abs(fit1(i)-fit2(i))
        if F<Fmin:
            Fmin=F
            mubest=i
            Tbest=0.5*(fit1(i)+fit2(i))
    return mubest,Tbest

def find_intersect_v1(mu_t1,mu_t2):             #fit in T mu space
    mu1=CubicSpline(mu_t1[0],mu_t1[1],bc_type=('not-a-knot',(2, 0.0)))
    mu2=CubicSpline(mu_t2[0],mu_t2[1],bc_type=('not-a-knot',(2, 0.0)))
    T=mu_t1[0][-1]-(mu_t1[1][-1]-mu_t2[1][-1])/(mu1(mu_t1[0][-1],1)-mu2(mu_t2[0][-1],1))
    mu=mu_t1[1][-1]+(T-mu_t1[0][-1])*mu1(mu_t1[0][-1],1)
    return T,mu

def plot_phase_boundary(xt):
    endslope=10
    w=np.ones(len(xt[0]))
    w[-1]=10
    if xt[0][-1]<xt[0][-2]:
        endslope=-10
    phb=UnivariateSpline(xt[1],xt[0],w=w)
    Tspace=np.linspace(xt[1][0],xt[1][-1],1000,endpoint=True)
    plt.plot(phb(Tspace),Tspace,'r')
    #plt.show()