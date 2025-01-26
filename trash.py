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

        signallist=np.array([2,0,0,0])
        print("shell mode, input whatever you want")
        with open("11.txt", "r") as file1:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt = eval(file1.read())
        plot_phase_diagram(xt,signallist,0.05,1,False,'g')
        print("check!!!!!!!!!!!!!!!!!!!!!!!!")
        print(signallist)
        with open("10.txt", "r") as file2:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt2 = eval(file2.read())
        plot_phase_diagram(xt2,signallist,0.05,1,False,'r')
        with open("09.txt", "r") as file3:
        # Use eval to parse the file content as a Python list with numpy arrays
            xt3 = eval(file3.read())
        plot_phase_diagram(xt3,signallist,0.05,1,False,'b')
        plt.xlim(0.1, 0.5)
        plt.ylim(1, 2.4)
        plt.xlabel("x")
        plt.ylabel("normalized T")
        plt.plot(0,0,color='g',label="r=1.1")
        plt.plot(0,0,color='r',label="r=1")
        plt.plot(0,0,color='b',label="r=0.9")
        plt.legend()
        plt.show()
    



def compute_ternary_cluster_energy(np.ndarray[DTYPE_t, ndim=4] prob, str lattice="BCC"):
    cdef DTYPE_t para
    if lattice == "BCC":
        para = -3.0
    else:
        para = 0.0

    cdef DTYPE_t abc000, abc001, abc010, abc011, abc100, abc101, abc110, abc111
    cdef DTYPE_t abd000, abd001, abd010, abd011, abd100, abd101, abd110, abd111
    cdef DTYPE_t acd000, acd001, acd010, acd011, acd100, acd101, acd110, acd111
    cdef DTYPE_t bcd000, bcd001, bcd010, bcd011, bcd100, bcd101, bcd110, bcd111

    # Compute sums
    abc000 = prob[0, 0, 0, 0] + prob[0, 0, 0, 1]
    abc001 = prob[0, 0, 1, 0] + prob[0, 0, 1, 1]
    abc010 = prob[0, 1, 0, 0] + prob[0, 1, 0, 1]
    abc011 = prob[0, 1, 1, 0] + prob[0, 1, 1, 1]
    abc100 = prob[1, 0, 0, 0] + prob[1, 0, 0, 1]
    abc101 = prob[1, 0, 1, 0] + prob[1, 0, 1, 1]
    abc110 = prob[1, 1, 0, 0] + prob[1, 1, 0, 1]
    abc111 = prob[1, 1, 1, 0] + prob[1, 1, 1, 1]

    abd000 = prob[0, 0, 0, 0] + prob[0, 0, 1, 0]
    abd001 = prob[0, 0, 0, 1] + prob[0, 0, 1, 1]
    abd010 = prob[0, 1, 0, 0] + prob[0, 1, 1, 0]
    abd011 = prob[0, 1, 0, 1] + prob[0, 1, 1, 1]
    abd100 = prob[1, 0, 0, 0] + prob[1, 0, 1, 0]
    abd101 = prob[1, 0, 0, 1] + prob[1, 0, 1, 1]
    abd110 = prob[1, 1, 0, 0] + prob[1, 1, 1, 0]
    abd111 = prob[1, 1, 0, 1] + prob[1, 1, 1, 1]

    acd000 = prob[0, 0, 0, 0] + prob[0, 1, 0, 0]
    acd001 = prob[0, 0, 0, 1] + prob[0, 1, 0, 1]
    acd010 = prob[0, 0, 1, 0] + prob[0, 1, 1, 0]
    acd011 = prob[0, 0, 1, 1] + prob[0, 1, 1, 1]
    acd100 = prob[1, 0, 0, 0] + prob[1, 1, 0, 0]
    acd101 = prob[1, 0, 0, 1] + prob[1, 1, 0, 1]
    acd110 = prob[1, 0, 1, 0] + prob[1, 1, 1, 0]
    acd111 = prob[1, 0, 1, 1] + prob[1, 1, 1, 1]

    bcd000 = prob[0, 0, 0, 0] + prob[1, 0, 0, 0]
    bcd001 = prob[0, 0, 0, 1] + prob[1, 0, 0, 1]
    bcd010 = prob[0, 0, 1, 0] + prob[1, 0, 1, 0]
    bcd011 = prob[0, 0, 1, 1] + prob[1, 0, 1, 1]
    bcd100 = prob[0, 1, 0, 0] + prob[1, 1, 0, 0]
    bcd101 = prob[0, 1, 0, 1] + prob[1, 1, 0, 1]
    bcd110 = prob[0, 1, 1, 0] + prob[1, 1, 1, 0]
    bcd111 = prob[0, 1, 1, 1] + prob[1, 1, 1, 1]

    # Compute the sum of x * log(x) for each term
    cdef DTYPE_t result = (
        abc000 * np.log(abc000) + abc001 * np.log(abc001) + abc010 * np.log(abc010) + abc011 * np.log(abc011) +
        abc100 * np.log(abc100) + abc101 * np.log(abc101) + abc110 * np.log(abc110) + abc111 * np.log(abc111) +
        abd000 * np.log(abd000) + abd001 * np.log(abd001) + abd010 * np.log(abd010) + abd011 * np.log(abd011) +
        abd100 * np.log(abd100) + abd101 * np.log(abd101) + abd110 * np.log(abd110) + abd111 * np.log(abd111) +
        acd000 * np.log(acd000) + acd001 * np.log(acd001) + acd010 * np.log(acd010) + acd011 * np.log(acd011) +
        acd100 * np.log(acd100) + acd101 * np.log(acd101) + acd110 * np.log(acd110) + acd111 * np.log(acd111) +
        bcd000 * np.log(bcd000) + bcd001 * np.log(bcd001) + bcd010 * np.log(bcd010) + bcd011 * np.log(bcd011) +
        bcd100 * np.log(bcd100) + bcd101 * np.log(bcd101) + bcd110 * np.log(bcd110) + bcd111 * np.log(bcd111)
    )

    return result * para

def compute_binary_cluster_energy(np.ndarray[DTYPE_t, ndim=4] prob, str lattice="BCC"):
    cdef np.ndarray[DTYPE_t, ndim=1] para
    if lattice == "BCC":
        para = np.array([1.5, 1], dtype=DTYPE)
    elif lattice == "FCC":
        para = np.array([-1, -1], dtype=DTYPE)
    else:
        raise ValueError("Invalid lattice type. Choose 'BCC' or 'FCC'.")

    cdef DTYPE_t ab00, ab01, ab10, ab11
    cdef DTYPE_t ac00, ac01, ac10, ac11
    cdef DTYPE_t ad00, ad01, ad10, ad11
    cdef DTYPE_t bc00, bc01, bc10, bc11
    cdef DTYPE_t bd00, bd01, bd10, bd11
    cdef DTYPE_t cd00, cd01, cd10, cd11

    # Compute sums
    ab00 = prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[0, 0, 1, 0] + prob[0, 0, 1, 1]
    ab01 = prob[0, 1, 0, 0] + prob[0, 1, 0, 1] + prob[0, 1, 1, 0] + prob[0, 1, 1, 1]
    ab10 = prob[1, 0, 0, 0] + prob[1, 0, 0, 1] + prob[1, 0, 1, 0] + prob[1, 0, 1, 1]
    ab11 = prob[1, 1, 0, 0] + prob[1, 1, 0, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1]

    ac00 = prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[0, 1, 0, 0] + prob[0, 1, 0, 1]
    ac01 = prob[0, 0, 1, 0] + prob[0, 0, 1, 1] + prob[0, 1, 1, 0] + prob[0, 1, 1, 1]
    ac10 = prob[1, 0, 0, 0] + prob[1, 0, 0, 1] + prob[1, 1, 0, 0] + prob[1, 1, 0, 1]
    ac11 = prob[1, 0, 1, 0] + prob[1, 0, 1, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1]

    ad00 = prob[0, 0, 0, 0] + prob[0, 0, 1, 0] + prob[0, 1, 0, 0] + prob[0, 1, 1, 0]
    ad01 = prob[0, 0, 0, 1] + prob[0, 0, 1, 1] + prob[0, 1, 0, 1] + prob[0, 1, 1, 1]
    ad10 = prob[1, 0, 0, 0] + prob[1, 0, 1, 0] + prob[1, 1, 0, 0] + prob[1, 1, 1, 0]
    ad11 = prob[1, 0, 0, 1] + prob[1, 0, 1, 1] + prob[1, 1, 0, 1] + prob[1, 1, 1, 1]

    bc00 = prob[0, 0, 0, 0] + prob[0, 0, 0, 1] + prob[1, 0, 0, 0] + prob[1, 0, 0, 1]
    bc01 = prob[0, 0, 1, 0] + prob[0, 0, 1, 1] + prob[1, 0, 1, 0] + prob[1, 0, 1, 1]
    bc10 = prob[0, 1, 0, 0] + prob[0, 1, 0, 1] + prob[1, 1, 0, 0] + prob[1, 1, 0, 1]
    bc11 = prob[0, 1, 1, 0] + prob[0, 1, 1, 1] + prob[1, 1, 1, 0] + prob[1, 1, 1, 1]

    bd00 = prob[0, 0, 0, 0] + prob[0, 0, 1, 0] + prob[1, 0, 0, 0] + prob[1, 0, 1, 0]
    bd01 = prob[0, 0, 0, 1] + prob[0, 0, 1, 1] + prob[1, 0, 0, 1] + prob[1, 0, 1, 1]
    bd10 = prob[0, 1, 0, 0] + prob[0, 1, 1, 0] + prob[1, 1, 0, 0] + prob[1, 1, 1, 0]
    bd11 = prob[0, 1, 0, 1] + prob[0, 1, 1, 1] + prob[1, 1, 0, 1] + prob[1, 1, 1, 1]

    cd00 = prob[0, 0, 0, 0] + prob[0, 1, 0, 0] + prob[1, 0, 0, 0] + prob[1, 1, 0, 0]
    cd01 = prob[0, 0, 0, 1] + prob[0, 1, 0, 1] + prob[1, 0, 0, 1] + prob[1, 1, 0, 1]
    cd10 = prob[0, 0, 1, 0] + prob[0, 1, 1, 0] + prob[1, 0, 1, 0] + prob[1, 1, 1, 0]
    cd11 = prob[0, 0, 1, 1] + prob[0, 1, 1, 1] + prob[1, 0, 1, 1] + prob[1, 1, 1, 1]

    # Compute results
    cdef DTYPE_t result1 = (
        ab00 * np.log(ab00) + ab01 * np.log(ab01) + ab10 * np.log(ab10) + ab11 * np.log(ab11) +
        cd00 * np.log(cd00) + cd01 * np.log(cd01) + cd10 * np.log(cd10) + cd11 * np.log(cd11)
    )

    cdef DTYPE_t result2 = (
        ac00 * np.log(ac00) + ac01 * np.log(ac01) + ac10 * np.log(ac10) + ac11 * np.log(ac11) +
        ad00 * np.log(ad00) + ad01 * np.log(ad01) + ad10 * np.log(ad10) + ad11 * np.log(ad11) +
        bc00 * np.log(bc00) + bc01 * np.log(bc01) + bc10 * np.log(bc10) + bc11 * np.log(bc11) +
        bd00 * np.log(bd00) + bd01 * np.log(bd01) + bd10 * np.log(bd10) + bd11 * np.log(bd11)
    )

    return para[0] * result1 + para[1] * result2

    if inputs.testmode:
        task="compute"
        composition=[0.38225,0.61775]
    
        print("code test mod")
        myCVM=CVM_loader("FCC",inputs,control_dict)
        #myCVM=FYLCVM(inputs.number_of_component,inputs.max_clustersize,inputs.basic_cluster_energy,inputs.vibration_parameter,inputs.local_energy_parameter,inputs.elastic_parameter,control_dict)                                #maybe considering using a class creator later
        if task=="point":
            print("compute the energy of given point")
            result=myCVM.optimize_free_energy(1,composition,"basinhopping")
            myCVM.output_optimization(result,1,composition)
        elif task=="trace":
            dictin={'range': [7.4,7.6], 'T': 1, 'phase': ['L10', 'L12'],'composition':[0.42,0.38]}
            phb=utility.phaseboundary(1,7.4,7.6,'L10', 'L12',0.42,0.38)
            #(newphb,muT,starting_pointnumber)=myCVM.trace_phase_boundary_v1(dictin)
            newphb=myCVM.trace_phase_boundary_v2(phb)
            print(phb.x1mat)
            print(phb.x2mat)
            print(phb.Tspace)
            #(newphb,muT,starting_pointnumber)=myCVM.trace_phase_boundary(dictin)
            #print(myCVM.starting_point_list)
            print(len(myCVM.node_list))
        elif task=="compute":
            #myCVM.compute_phase_diagram(control_dict['Tstart'],0,25,inputs.phasediagram_name)
            myCVM.compute_phase_diagram_v2(control_dict['Tstart'],0,10,inputs.phasediagram_name)
        elif task=="scatter":
            #myCVM.scan_phase_diagram(control_dict['Tstart'],-75,75)
            myCVM.scan_phase_diagram(1,0,10)
        elif task=="scan":
            myCVM.scan_phase_boundary(1,0,10)
        elif task=="find":
            myCVM.find_phase_boundary(500,900,[0.5,0.5],"brute")
        elif task=="brute":
            inputlist=[(10,-20,"L10"),(10,-10,"L10"),(10,0,"L12"),(10,10,"L12")]
            for inputvalue in inputlist:
                myCVM.plot_potential_surface(inputvalue[0],inputvalue[1],inputvalue[2])
        elif task=="area":
            #mu is between 24.23 and 24.25 composition is between 0.408 and 0.389 at T=500.00
            Trange=np.array([490,520])
            murange=np.array([23,45])
            myCVM.plot_point_region(Trange,murange,5,0.1)
        elif task=="test":
            print("just check if input is correct")
            print(inputs.phasediagram_name)
        print("--- %s seconds ---" % (time.time() - start_time))