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

            testinput=[np.array([[8.35311071e-01, 8.34217126e-01, 8.33236638e-01, 8.32144872e-01,
        8.31055195e-01, 8.29966810e-01, 8.28880412e-01, 8.27684699e-01,
        8.26490753e-01, 8.25300744e-01, 8.24110460e-01, 8.22922697e-01,
        8.21628570e-01, 8.20336231e-01, 8.19047684e-01, 8.17652946e-01,
        8.16260950e-01, 8.14872293e-01, 8.13379752e-01, 8.11889819e-01,
        8.10403851e-01, 8.08813723e-01, 8.07227628e-01, 8.05540316e-01,
        8.03856291e-01, 8.02070982e-01, 8.00289607e-01, 7.98409745e-01,
        7.96429189e-01, 7.94349985e-01, 7.92275419e-01, 7.89999165e-01,
        7.87730740e-01, 7.85261121e-01, 7.82695786e-01, 7.79932696e-01,
        7.76974435e-01, 7.73718015e-01, 7.70168619e-01, 7.66223053e-01,
        7.61687466e-01, 7.55961134e-01, 7.48055799e-01],
       [7.87384561e-01, 7.86697673e-01, 7.86210212e-01, 7.85576828e-01,
        7.84965414e-01, 7.84384053e-01, 7.83822320e-01, 7.83137287e-01,
        7.82472311e-01, 7.81838924e-01, 7.81231471e-01, 7.80647497e-01,
        7.79953196e-01, 7.79284852e-01, 7.78643454e-01, 7.77902534e-01,
        7.77191972e-01, 7.76505745e-01, 7.75734836e-01, 7.74996229e-01,
        7.74286411e-01, 7.73486758e-01, 7.72726079e-01, 7.71892058e-01,
        7.71085856e-01, 7.70214109e-01, 7.69371334e-01, 7.68471274e-01,
        7.67513437e-01, 7.66507650e-01, 7.65527617e-01, 7.64436422e-01,
        7.63359833e-01, 7.62175567e-01, 7.60944769e-01, 7.59607218e-01,
        7.58162240e-01, 7.56553258e-01, 7.54778115e-01, 7.52764298e-01,
        7.50397038e-01, 7.47298437e-01, 7.42787530e-01],
       [5.50000000e+02, 5.55000000e+02, 5.60000000e+02, 5.65000000e+02,
        5.70000000e+02, 5.75000000e+02, 5.80000000e+02, 5.85000000e+02,
        5.90000000e+02, 5.95000000e+02, 6.00000000e+02, 6.05000000e+02,
        6.10000000e+02, 6.15000000e+02, 6.20000000e+02, 6.25000000e+02,
        6.30000000e+02, 6.35000000e+02, 6.40000000e+02, 6.45000000e+02,
        6.50000000e+02, 6.55000000e+02, 6.60000000e+02, 6.65000000e+02,
        6.70000000e+02, 6.75000000e+02, 6.80000000e+02, 6.85000000e+02,
        6.90000000e+02, 6.95000000e+02, 7.00000000e+02, 7.05000000e+02,
        7.10000000e+02, 7.15000000e+02, 7.20000000e+02, 7.25000000e+02,
        7.30000000e+02, 7.35000000e+02, 7.40000000e+02, 7.45000000e+02,
        7.50000000e+02, 7.55000000e+02, 7.60000000e+02]]), np.array([[6.04130640e-01, 6.02834946e-01, 6.01547495e-01, 6.00153219e-01,
        5.98652148e-01, 5.97052100e-01, 5.95349526e-01, 5.93549259e-01,
        5.91648408e-01, 5.89541855e-01],
       [5.87149177e-01, 5.86480492e-01, 5.85852091e-01, 5.85066530e-01,
        5.84117150e-01, 5.82998554e-01, 5.81728337e-01, 5.80276669e-01,
        5.78661021e-01, 5.76652794e-01],
       [5.50000000e+02, 5.55000000e+02, 5.60000000e+02, 5.65000000e+02,
        5.70000000e+02, 5.75000000e+02, 5.80000000e+02, 5.85000000e+02,
        5.90000000e+02, 5.95000000e+02]]), np.array([[4.29352418e-01, 4.31445971e-01, 4.33452983e-01, 4.35624736e-01,
        4.37706224e-01, 4.39705068e-01, 4.41624890e-01, 4.43663659e-01,
        4.45617108e-01, 4.47683690e-01, 4.49464897e-01, 4.51527910e-01,
        4.53661926e-01, 4.56941511e-01, 4.60126886e-01, 4.63339926e-01,
        4.66534024e-01, 4.69805220e-01, 4.73092063e-01, 4.76357482e-01,
        4.79730492e-01, 4.82850934e-01, 4.86857370e-01, 4.89645308e-01,
        4.94451239e-01, 4.96190752e-01, 5.01680149e-01, 5.01563307e-01,
        5.01709529e-01],
       [4.04286718e-01, 4.05978205e-01, 4.07665713e-01, 4.09464599e-01,
        4.11260744e-01, 4.13050885e-01, 4.14837365e-01, 4.16733546e-01,
        4.18626562e-01, 4.20626093e-01, 4.22508458e-01, 4.24612208e-01,
        4.26823548e-01, 4.29925561e-01, 4.33127261e-01, 4.36538870e-01,
        4.40158348e-01, 4.44093491e-01, 4.48340413e-01, 4.52895501e-01,
        4.57972346e-01, 4.63133079e-01, 4.70201279e-01, 4.75734193e-01,
        4.85593780e-01, 4.89476189e-01, 5.01025151e-01, 5.00676148e-01,
        5.00853388e-01],
       [5.50000000e+02, 5.55000000e+02, 5.60000000e+02, 5.65000000e+02,
        5.70000000e+02, 5.75000000e+02, 5.80000000e+02, 5.85000000e+02,
        5.90000000e+02, 5.95000000e+02, 6.00000000e+02, 6.05000000e+02,
        6.10000000e+02, 6.15000000e+02, 6.20000000e+02, 6.25000000e+02,
        6.30000000e+02, 6.35000000e+02, 6.40000000e+02, 6.45000000e+02,
        6.50000000e+02, 6.55000000e+02, 6.60000000e+02, 6.65000000e+02,
        6.70000000e+02, 6.75000000e+02, 6.80000000e+02, 6.85000000e+02,
        6.90000000e+02]]), np.array([[5.88081448e-01, 5.95477038e-01, 6.01642058e-01, 6.07327825e-01,
        6.12396430e-01, 6.17055761e-01, 6.21532290e-01, 6.25686053e-01,
        6.29757878e-01, 6.33740684e-01, 6.37491033e-01, 6.41133241e-01,
        6.44789610e-01, 6.48461166e-01, 6.51986247e-01, 6.55519429e-01,
        6.59034721e-01, 6.62532613e-01, 6.66000757e-01, 6.69438033e-01,
        6.72964916e-01, 6.76434199e-01, 6.79978461e-01, 6.83575104e-01,
        6.87207277e-01, 6.90988891e-01, 6.94769566e-01, 6.98644364e-01,
        7.02806120e-01, 7.07182944e-01, 7.11914400e-01, 7.17252348e-01,
        7.24050753e-01],
       [5.87444449e-01, 5.94592094e-01, 6.00496354e-01, 6.05876289e-01,
        6.10632472e-01, 6.14969839e-01, 6.19093590e-01, 6.22903294e-01,
        6.26603617e-01, 6.30194142e-01, 6.33575748e-01, 6.36846811e-01,
        6.40114099e-01, 6.43373828e-01, 6.46526644e-01, 6.49673319e-01,
        6.52814287e-01, 6.55950129e-01, 6.59081018e-01, 6.62205420e-01,
        6.65426841e-01, 6.68641747e-01, 6.71952837e-01, 6.75359543e-01,
        6.78860928e-01, 6.82555001e-01, 6.86346132e-01, 6.90330221e-01,
        6.94708398e-01, 6.99478941e-01, 7.04839417e-01, 7.11188992e-01,
        7.19817720e-01],
       [6.00000000e+02, 6.05000000e+02, 6.10000000e+02, 6.15000000e+02,
        6.20000000e+02, 6.25000000e+02, 6.30000000e+02, 6.35000000e+02,
        6.40000000e+02, 6.45000000e+02, 6.50000000e+02, 6.55000000e+02,
        6.60000000e+02, 6.65000000e+02, 6.70000000e+02, 6.75000000e+02,
        6.80000000e+02, 6.85000000e+02, 6.90000000e+02, 6.95000000e+02,
        7.00000000e+02, 7.05000000e+02, 7.10000000e+02, 7.15000000e+02,
        7.20000000e+02, 7.25000000e+02, 7.30000000e+02, 7.35000000e+02,
        7.40000000e+02, 7.45000000e+02, 7.50000000e+02, 7.55000000e+02,
        7.60000000e+02]]), np.array([[5.86726059e-01, 5.84320554e-01, 5.81917907e-01, 5.79414033e-01,
        5.76911149e-01, 5.74407871e-01, 5.71906919e-01, 5.69198981e-01,
        5.66185822e-01, 5.61117900e-01, 5.56047308e-01, 5.50870803e-01,
        5.44866087e-01, 5.37098433e-01, 5.28381387e-01, 5.25759011e-01,
        5.14297742e-01, 5.03087129e-01, 5.02320981e-01],
       [5.74252175e-01, 5.71067449e-01, 5.67967560e-01, 5.64769224e-01,
        5.61683749e-01, 5.58716364e-01, 5.55864674e-01, 5.52807740e-01,
        5.49425811e-01, 5.43372926e-01, 5.37970402e-01, 5.33020703e-01,
        5.27801575e-01, 5.21776117e-01, 5.15988284e-01, 5.14713443e-01,
        5.08306705e-01, 5.02668544e-01, 5.02356514e-01],
       [6.00000000e+02, 6.05000000e+02, 6.10000000e+02, 6.15000000e+02,
        6.20000000e+02, 6.25000000e+02, 6.30000000e+02, 6.35000000e+02,
        6.40000000e+02, 6.45000000e+02, 6.50000000e+02, 6.55000000e+02,
        6.60000000e+02, 6.65000000e+02, 6.70000000e+02, 6.75000000e+02,
        6.80000000e+02, 6.85000000e+02, 6.90000000e+02]])]