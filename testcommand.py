#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

x1=np.array([1,2,3,4,5,6])
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
print(x)