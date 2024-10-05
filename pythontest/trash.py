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