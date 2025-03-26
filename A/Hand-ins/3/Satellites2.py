import numpy as np
import matplotlib.pyplot as plt
import copy 

def n(x,A,Nsat,a,b,c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def N(x, A, Nsat, a, b, c):
    return n(x, A, Nsat, a, b, c)*4*np.pi*x**2

a = 2.4
b = 0.25
c = 1.6
Nsat = 100
A = 256/(5*np.pi**(3/2))

def crout(A):
    m, n = A.shape
    for i in range(m):
        for j in range(n):
            if i >= j: # if i > j, add the alpha[i, j] element
                A[i, j] = A[i, j] - np.sum(A[i, :j]*A[:j, j])
            else: # else, edit matrix beta
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

def f(x):
    return (((x+10)*x + 10)*x - 40)*x + 40

phi = (1+np.sqrt(5)) * 0.5

# Construct the Vandermonde matrix.
def solve_system(func, a, b, c, *args):
    A = np.zeros([3, 3])
    x = np.array([a, b, c])
    y = func(x, *args)

    m, n = A.shape
    for i in range(m):
        for j in range(n):
            A[i, j] = x[i]**j
    A = np.fliplr(A)
    coeffs, _ = get_coeffs(y, A)

    return coeffs


def find_bracket(func, brac, *args):
    a = brac[0]
    b = brac[1]
  

    if func(a, *args) < func(b, *args):
        print('swapped a and b')
        a = brac[1]
        b = brac[0]
    
    c = b + (b-a)*phi
    while func(c, *args) < func(b, *args):
        # let's fit a parabola. 
        p, q, r = solve_system(func, a, b, c, *args)
        # plt.plot(x, p*x*x + q*x + r)
        d = -q/(2*p)
        if func(d, *args) < func(c, *args):
            return [b, d, c]
        elif func(d, *args) > func(b, *args):
            return [a, b, d]
        else: 
            d = c + (c-b)*phi
        if abs(d-b) > 100*abs(c-b):
            d = c + (c-b)*phi
        a = b
        b = c
        c = d 
    return [a, b, c]


def golden_section(func, brac, *args, acc = 1e-5):
    a, b, c = brac[0], brac[1], brac[2]
    w = 2 - phi 

    while abs(c-a) > acc:       
        if abs(c-b) > abs(b-a):
            d = b + (c - b)*w
            if func(d, *args) < func(b, *args):
                a = b
                b = d
            else:
                c = d
        else:
            d = b + (a - b)*w
            if func(d, *args) < func(b, *args):
                c = b
                b = d
            else:
                a = d
    if func(d, *args) < func(b, *args):
        return d
    else: 
        return b

xa = np.linspace(0.0001, 5, 1000)
Na = N(xa, A, Nsat, a, b, c)

def negative_N(xa, A, Nsat, a, b, c):
    return -1 * N(xa, A, Nsat, a, b, c)

brac = [1, 5]
brac = find_bracket(negative_N, brac, A, Nsat, a, b, c)

max_x = golden_section(negative_N, brac, A, Nsat, a, b, c)

print('x value at maximum:', max_x, '\nN(x) value at maximum:', N(max_x, A, Nsat, a, b, c))


def readfile(filename):
    f = open(filename, 'r')
    data = f.readlines()[3:] #Skip first 3 lines 
    nhalo = int(data[0]) #number of halos
    radius = []
    
    for line in data[1:]:
        if line[:-1]!='#':
            radius.append(float(line.split()[0]))
    
    radius = np.array(radius, dtype=float)    
    f.close()
    return radius, nhalo #Return the virial radius for all the satellites in the file, and the number of halos

#Call this function as: 
#radius, nhalo = readfile('satgals_m15.txt')


# Plot of binned data with the best fit (question 1b.4 and 1c)
# As always, feel free to replace by your own plotting routines if you want
xmin, xmax = 1e-4, 5. # replace by your choices
n_bins = 100 # replace by your binning
edges = np.exp(np.linspace(np.log(xmin), np.log(xmax), n_bins+1))

fig1b, ax = plt.subplots(3,2,figsize=(6.4,8.0))
for i in range(5):
    Nsat = 100 # replace by actual appropriate number for mass bin i
    x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for mass bin i
    Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for mass bin i integrated per radial bin
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
    row=i//2
    col=i%2
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].step(edges[:-1], Ntilda, where='post', label='best-fit profile')
    ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1b.png', dpi=600)

# Plot 1c (same code as above)
fig1c, ax = plt.subplots(3,2,figsize=(6.4,8.0))
for i in range(5):
    Nsat = 100 # replace by actual appropriate number for mass bin i
    x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for mass bin i
    Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for mass bin i integrated per radial bin
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
    row=i//2
    col=i%2
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].step(edges[:-1], Ntilda, where='post', label='best-fit profile')
    ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1c.png', dpi=600)


# BONUS: Monte Carlo resampled fits (1e)
num_samples = 10 #replace by how many resamplings you can draw/fit in reasonable time
fig1e, ax = plt.subplots()
Nsat = 100 # replace by actual appropriate number for mass bin i
x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for chosen mass bin
binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
ax.step(edges[:-1], binned_data, where='post', label='binned data')
Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
ax.step(edges[:-1], Ntilda, where='post', label='best-fit profiles', color="C1")
for i in range(num_samples):
    Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
    ax.step(edges[:-1], Ntilda, where='post', color="C1")
# Also plot the mean or median fitted profile here
ax.set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{...}} M_{{\\odot}}/h$")
plt.legend()
plt.savefig('my_solution_1e_chisq.png', dpi=600)

num_samples = 10 #replace by how many resamplings you can draw/fit in reasonable time
fig1e, ax = plt.subplots()
Nsat = 100 # replace by actual appropriate number for mass bin i
x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for chosen mass bin
binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
ax.step(edges[:-1], binned_data, where='post', label='binned data')
Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
ax.step(edges[:-1], Ntilda, where='post', label='best-fit profiles', color="C2")
for i in range(num_samples):
    Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
    ax.step(edges[:-1], Ntilda, where='post', color="C2")
# Also plot the mean or median fitted profile here
ax.set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{...}} M_{{\\odot}}/h$")
plt.legend()
plt.savefig('my_solution_1e_poisson.png', dpi=600)