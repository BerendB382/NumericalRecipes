# Main code body and plotting code adapted from Marcel van Daalen, at
# https://home.strw.leidenuniv.nl/~daalen/Handin_files/vandermonde.py

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import copy
from time import perf_counter

os.chdir('A/Hand-ins/1')

def crout(A):
    m, n = A.shape
    for i in range(m):
        for j in range(n):
            if i >= j:
                A[i, j] = A[i, j] - np.sum(A[i, :j]*A[:j, j])
            else:
                A[i, j] = (A[i, j] - np.sum(A[i, :i]*A[:i, j]))/A[i, i]
        
    return A

def forward_sub2(A, b):
    m, _ = A.shape
    y = copy.deepcopy(b)
    for k in range(m):
        y[k] = (b[k] - np.sum(A[k, :k]*y[:k]))/A[k, k]
    return y

def backward_sub2(A, b):
    _, n = A.shape
    x = np.zeros_like(b)
    for i in range(n-1, -1, -1):
        x[i] = b[i] - np.dot(A[i, i+1:], x[i+1:])
    return x

def get_coeffs(y, A):
    A_decomp = crout(A)
    fw_sub = forward_sub2(A_decomp, y)
    coeffs = backward_sub2(A_decomp, fw_sub)
    
    return coeffs, A_decomp

def polynomial(x, coeffs):
    y = np.zeros_like(x)
    for i in range(len(x)):
        for j in range(len(coeffs)):
            y[i] += coeffs[j] * x[i]**j
    return y

# Let's get started 
data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)

x=data[:,0]
y=data[:,1]

# Construct the matrix.
A = np.zeros([20, 20])
m, n = A.shape
for i in range(m):
    for j in range(n):
        A[i, j] = x[i]**j

A_orig = copy.deepcopy(A)
# now set up the polynomial function!
B = np.zeros([5, 5])
for i in range(5):
    for j in range(5):
        B[i, j] = x[i]**j

coeffs, crout_LU = get_coeffs(y, A)
xx=np.linspace(x[0],x[-1],1001) #x values to interpolate at
yya = polynomial(xx, coeffs)

# Check crout 
# A_decomp = crout(A)
# U = np.triu(A_decomp, k=0)
# L = np.tril(A_decomp, k = -1) + np.identity(len(A_decomp))

# print(A_orig - (L @ U))

# Now we'd like to implement neville's algorithm.

def neville(datax, datay, x_interp):
    y_interp = []
    error_est = []
    for x in x_interp:
        y = datay.copy()
        M = len(datax)

        for k in range(1, M):
            for i in range(M-k):
                y[i] = ((x - datax[i+k])*y[i] + (datax[i] - x)*y[i+1]) / (datax[i] - datax[i+k])
        y_interp.append(y[0])
        error_est.append(np.abs(y[0] - y[1]))
    return y_interp, error_est 
    
# Let's implement the error canceling algorithm.

def error_cancel(A_orig, crout_LU, y, coeffs, iterations):
    for i in range(iterations):
        v = A_orig @ coeffs - y
        fw_result = forward_sub2(crout_LU, v)
        coeff_corr = backward_sub2(crout_LU, fw_result)
        coeffs = coeffs - coeff_corr
    return coeffs

#The plotting code below assumes you've given the interpolated
#values for 2a suffix "a", those for 2b "b", and those for 2c 
#"c1" and "c10" – feel free to change these
#Note that you interpolate to xx for the top panel but to
#x for the bottom panel (since we can only compare the accuracy
#to the given points, not in between them)
ya= polynomial(x, coeffs)
yyb, yyb_err= neville(x, y, xx)
yb, yb_err= neville(x, y, x)

c1_coeffs = error_cancel(A_orig, crout_LU, y, coeffs, iterations=1)
yyc1= polynomial(xx, c1_coeffs) 
yc1= polynomial(x, c1_coeffs)
c10_coeffs = error_cancel(A_orig, crout_LU, y, coeffs, iterations=10)

yyc10= polynomial(xx, c10_coeffs)
yc10= polynomial(x, c10_coeffs)

#Don't forget to output the coefficients you find with your LU routine
print('Coefficients found with LU method (a):\n', coeffs)

print('Coefficients found with 1 iteration of LU correction (c):\n', c1_coeffs)

print('Coefficients found with 10 iterations of LU correction (c):\n', c10_coeffs)
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