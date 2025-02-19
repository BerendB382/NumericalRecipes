# Main code body and plotting code adapted from Marcel van Daalen, at
# https://home.strw.leidenuniv.nl/~daalen/Handin_files/vandermonde.py

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

def crout(A):
    '''Overwrites a matrix by its LU composition using Crout's algorithm.'''
    assert A.shape == A.T.shape
    alpha = np.identity(len(A))
    for j in range(len(A)):
        for i in range(len(A)):
            if i <= j:
                if i >= 0:
                    sumvar_a = np.arange(0, i, 1)
                    sum_a = np.sum(alpha[i, sumvar_a]*A[sumvar_a, j])
                    A[i, j] -= sum_a
                
            else:
                if j > 0:
                    sumvar_b = np.arange(0, j, 1)
                else:
                    sumvar_b = np.array([int(0)]) 
                sum_b = np.sum(alpha[i, sumvar_b]*A[sumvar_b, j])
                A[i, j] = 1 / A[j, j] * (A[i, j] - sum_b)
                if i != j:
                    alpha[i, j] = A[i, j]
        
    return A



data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)
x=data[:,0]
y=data[:,1]
xx=np.linspace(x[0],x[-1],1001) #x values to interpolate at

#Insert your own code to calculate interpolated y values here!
#The plotting code below assumes you've given the interpolated
#values for 2a suffix "a", those for 2b "b", and those for 2c 
#"c1" and "c10" – feel free to change these
#Note that you interpolate to xx for the top panel but to
#x for the bottom panel (since we can only compare the accuracy
#to the given points, not in between them)
yya=np.zeros(1001,dtype=np.float64) #replace!
ya=np.zeros(len(x),dtype=np.float64) #replace!
yyb=yya #replace!
yb=ya #replace!
yyc1=yya #replace!
yc1=ya #replace!
yyc10=yya #replace!
yc10=ya #replace!

#Don't forget to output the coefficients you find with your LU routine

#Plot of points with absolute difference shown on a log scale (question 2a)
fig=plt.figure()
gs=fig.add_gridspec(2,hspace=0,height_ratios=[2.0,1.0])
axs=gs.subplots(sharex=True,sharey=False)
axs[0].plot(x,y,marker='o',linewidth=0)
plt.xlim(-1,101)
axs[0].set_ylim(-400,400)
axs[0].set_ylabel('$y$')
axs[1].set_ylim(1e-16,1e1)
axs[1].set_ylabel('$|y-y_i|$')
axs[1].set_xlabel('$x$')
axs[1].set_yscale('log')
line,=axs[0].plot(xx,yya,color='orange')
line.set_label('Via LU decomposition')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-ya),color='orange')
plt.savefig('my_vandermonde_sol_2a.png',dpi=600)

#For questions 2b and 2c, add this block
line,=axs[0].plot(xx,yyb,linestyle='dashed',color='green')
line.set_label('Via Neville\'s algorithm')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-yb),linestyle='dashed',color='green')
plt.savefig('my_vandermonde_sol_2b.png',dpi=600)

#For question 2c, add this block too
line,=axs[0].plot(xx,yyc1,linestyle='dotted',color='red')
line.set_label('LU with 1 iteration')
axs[1].plot(x,abs(y-yc1),linestyle='dotted',color='red')
line,=axs[0].plot(xx,yyc10,linestyle='dashdot',color='purple')
line.set_label('LU with 10 iterations')
axs[1].plot(x,abs(y-yc10),linestyle='dashdot',color='purple')
axs[0].legend(frameon=False,loc="lower left")
plt.savefig('my_vandermonde_sol_2c.png',dpi=600)

#Don't forget to caption your figures to describe them/
#mention what conclusions you draw from them!