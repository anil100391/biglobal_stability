import matplotlib.pyplot as plt
import numpy as np
from math import cos,pi
def plot_file(a,nx,ny):
    ythe = []
    ynum = []
    x = np.zeros(21)
    for i in range(0,nx+1):
        x[i] = cos(i*pi/nx)
    y = x
    for line in open(a):
        columns = line.split('\t')
        ythe.append(float(columns[0]))
        ynum.append(float(columns[1]))

    plt.plot(y,ythe[0:21],'b',label='Theoretical')
    plt.xlabel('y',fontsize=10,family='sans-serif')
    plt.ylabel('f(0,y)',fontsize=10,family='serif')
    plt.grid(True)
    plt.plot(y,ynum[0:21],'ko',label='Numerical')
    plt.legend(loc=1,family='sans-serif')
    plt.savefig('plt.png')

plot_file('data.txt',20,20)
    
