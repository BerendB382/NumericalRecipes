import numpy as np
import matplotlib.pyplot as plt
import copy 
import sys
import os
from pathlib import Path

def n(x,A,Nsat,a,b,c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def N(x, A, Nsat, a, b, c):
    return n(x, A, Nsat, a, b, c)*4*np.pi*x**2

def n_deriv_a(x,A,Nsat,a,b,c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c) * np.log((x/b))*4*np.pi*x**2

def n_deriv_b(x,A,Nsat,a,b,c):
    return (A*Nsat*(c*(x/b)**(c) - a + 3)*(x/b)**(a)*b**2 * np.exp(-(x/b)**c)) * x**(-3)  *4*np.pi*x**2

def n_deriv_c(x,A,Nsat,a,b,c):
    return -A*Nsat*((x/b)**(c+a-3))*np.log(x/b)*np.exp(-(x/b)**c)*4*np.pi*x**2

def n_deriv_A(x,A,Nsat,a,b,c):
    return Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)*4*np.pi*x**2

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

phi = (1+np.sqrt(5)) * 0.5

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

### 1B.

def romberg_interval(func, a, b, order, *args):
    '''Calculates the integral of a function from a to b, using 
    some order of Richardson extrapolation'''
    r = np.zeros(order)
    
    # calculate r0
    h = (b - a)
    r[0] = 0.5 * h * (func(a, *args) - func(b, *args))
    Np = 1
    for i in np.arange(1, order-1): # calculate approximations
        r[i] = 0
        delta = h
        h = 0.5 * h
        x = a + h
        for n in range(Np):
           r[i] += func(x, *args)
           x += delta
        r[i] = 0.5*(r[i-1] + delta*r[i])
        Np = 2*Np
    Np = 1
    for i in np.arange(1, order-1): # combine approximations. r[0] holds the best one
        Np = 4*Np 
        for j in np.arange(0, order-i):
            r[j] = (Np * r[j+1] - r[j]) / (Np - 1)

    return r

def romberg(func, x, order, *args):
    '''with x a linspace. calculates the integral over a range given by x.
    len(x) is the number of function evaluations we use to integrate.'''
    result = 0
    for i in range(len(x)-1):
        a = x[i]
        b  = x[i+1]
        result += romberg_interval(func, a, b, order, *args)[0]
    return result

def readfile(filename):
    os.chdir('/Users/bjhnieuwhof/Google Drive/Universiteit Leiden/Master Astronomy/Numerical Recipes/A/Hand-ins/3')
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

def get_A(x, Nsat, a, b, c):
    integ_result = romberg(N, x, 10, 1, Nsat, a, b, c)
    return Nsat / integ_result

def bin_data(radius, nhalo, xmin, xmax):
    nhaloTh = 1 / nhalo
    edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 21)
    binned_data, edges = np.histogram(radius, bins=edges)
    binwidths = np.diff(edges)
    N_per_halo = binned_data * nhaloTh / binwidths # normalize by dividing by the bin widths 

    # average Nsat
    Nsat = len(radius) / nhalo
    return N_per_halo, edges, Nsat

def z(x_data, y_data, sigma, model, *args):
    sigma_log = sigma / (y_data * np.log(10))
    return (np.log10(y_data) - np.log10(model(x_data, *args))) / sigma_log

def chi_squared(x_data, y_data, sigma, model, *args):
    return np.sum(z(x_data, y_data, sigma, model, *args)**2)

def chi(x_data, y_data, sigma, model, deriv, *args):
    sigma_log = sigma / (y_data * np.log(10))
    return np.sum(2*deriv(x_data, *args)*(z(x_data, y_data, sigma, model, *args)) / sigma_log)

def calculate_alpha(x_data, p, sigma, lamb):
    n_param = len(p)-2
    alpha = np.zeros((n_param, n_param))
    funclist = [n_deriv_a, n_deriv_b, n_deriv_c]
    for k in range(n_param):
        for j in range(n_param):
            alpha[k, j] = np.sum(funclist[k](x_data, *p)*funclist[j](x_data, *p) / sigma**2)
    for i in range(n_param):
        alpha[i, i] *= (1 + lamb)
    return alpha

def calculate_beta(x_data, y_data, p, sigma, func, weight):
    n_param = len(p)-2
    beta = np.zeros(n_param)
    xderivlist = [n_deriv_a, n_deriv_b, n_deriv_c]
    for k in range(n_param):
        beta[k] = 0.5*weight(x_data, y_data, sigma, func, xderivlist[k], *p) # * -1
    return beta

def levenberg_marquardt(x_data, y_data, initial_guess, sigma, func, metric_func, weight_func, lamb=1e-3, w = 10, acc = 0.1):
    # TODO: Calculate the metric in log space.
    metric = metric_func(x_data, y_data, sigma, func, *initial_guess)
    lamb = lamb
    j = 0
    p = initial_guess
    while True:
        alpha = calculate_alpha(x_data, p, sigma, lamb)
        beta = calculate_beta(x_data, y_data, p, sigma, func, weight_func)
        prev_p = np.copy(p)
        dp = np.zeros_like(prev_p) 
        dp[2:], _ = get_coeffs(beta, alpha) 
        new_p = prev_p + dp
        # recalculate A based on new a,b,c
        new_p[0] = get_A(x_data, Nsat, *new_p[2:])

        new_metric = metric_func(x_data, y_data, sigma, func, *new_p)

        delta_metric = abs(new_metric - metric)

        if new_metric >= metric:
            lamb *= w 
            metric = metric
            p = p
        else: 
            metric = new_metric
            p = new_p
            lamb /= w
        j += 1
        if j > 1000:
            print('max iterations reached')
            return p, j
        if abs(delta_metric/len(initial_guess)) < acc:
            print(f'accuracy threshold reached in {j} iterations')
            return p, j

def fit_data_set(filename, initial_guess, acc=0.001):
    xmin, xmax = 1e-4, 4
    radius, nhalo = readfile(filename)
    N_per_halo, edges, Nsat = bin_data(radius, nhalo, xmin, xmax)
    bin_centers = (edges[:-1] + edges[1:]) * 0.5
    
    plt.stairs(N_per_halo, edges, label='binned_data')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-4, 1e-1)
    plt.xlabel('bins (log)')
    plt.ylabel(f'Counts per halo (log)')
    plt.title(f'Loglog plot for {filename}, Nsat = {Nsat}')
    
    # create the input parameters.
    A = get_A(bin_centers, Nsat, *initial_guess)
    params_i = [A, Nsat, *initial_guess]

    # x data will be 
    p, j = levenberg_marquardt(bin_centers, N_per_halo, params_i, np.sqrt(N_per_halo), N, chi_squared, chi, acc=acc)
    print(p)
    test_x = 10**np.linspace(np.log10(xmin), np.log10(xmax), 100)
    plt.plot(test_x, N(test_x, *params_i), label='initial guess')
    plt.plot(test_x, N(test_x, *p), label='fit')
    plt.legend()
    plt.show()
    chi2_min = chi_squared(bin_centers, N_per_halo, np.sqrt(N_per_halo), n, *p)
    return Nsat, p, chi2_min, bin_centers, edges

#Call this function as: 
fit_data_set('satgals_m11.txt', [2.4, 0.25, 1.6], acc=1e-40)
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