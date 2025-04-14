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
brac = find_bracket(negative_N, brac, A, Nsat, a, b, c) # use negative N, as we're minimizing.

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
    edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 51)
    binned_data, edges = np.histogram(radius, bins=edges)
    binwidths = np.diff(edges)
    N_per_halo = binned_data * nhaloTh / binwidths # normalize by dividing by the bin widths and the number of halos

    # average Nsat
    Nsat = len(radius) / nhalo
    return N_per_halo, edges, Nsat

def edges_to_bins(edges):
    edges = np.asarray(edges) # function to make edges usable for romberg_interval
    return np.column_stack((edges[:-1], edges[1:]))

def N_integrated(edges, A, Nsat, a, b, c):
    result = np.zeros_like(edges[:,0])
    binwidths = np.diff(edges)
    inverse_binwidths = 1 / binwidths[:, 0]
    for i in range(len(result)):
        xmin, xmax = edges[i]
        result[i] = romberg_interval(N, xmin, xmax, 5, A, Nsat, a, b, c)[0] * inverse_binwidths[i] #normalize by the bin widths also
    return result 

def z(x_data, y_data, edges, sigma, model, *args):
    return (y_data - model(edges, *args)) / sigma

def chi_squared(x_data, y_data, edges, sigma, model, *args):
    return np.sum(0.5*z(x_data, y_data, edges, sigma, model, *args)**2)

def chi(x_data, y_data, edges, sigma, model, deriv, *args):
    return np.sum(deriv(x_data, *args)*(z(x_data, y_data, edges, sigma, model, *args)) / sigma)

def calculate_alpha(x_data, y_data, p, sigma, lamb):
    n_param = len(p)-2
    alpha = np.zeros((n_param, n_param))
    funclist = [n_deriv_a, n_deriv_b, n_deriv_c]
    for k in range(n_param):
        for j in range(n_param):
            alpha[k, j] = np.sum(funclist[k](x_data, *p)*funclist[j](x_data, *p) / sigma**2)
    for i in range(n_param):
        alpha[i, i] *= (1 + lamb)
    return alpha

def calculate_beta(x_data, y_data, edges, p, sigma, func, weight):
    n_param = len(p)-2
    beta = np.zeros(n_param)
    xderivlist = [n_deriv_a, n_deriv_b, n_deriv_c]
    for k in range(n_param):
        beta[k] = 0.5*weight(x_data, y_data, edges, sigma, func, xderivlist[k], *p)
    return beta

def levenberg_marquardt(x_data, y_data, edges, initial_guess, sigma, func, metric_func, weight_func, lamb=0.001, w = 10, acc = 0.1):
    bins = edges_to_bins(edges)
    metric = metric_func(x_data, y_data, bins, sigma, func, *initial_guess)
    lamb = lamb
    j = 0
    p = initial_guess
    while True:
        alpha = calculate_alpha(x_data, y_data, p, sigma, lamb)
        beta = calculate_beta(x_data, y_data, bins, p, sigma, func, weight_func)
        prev_p = np.copy(p)
        dp = np.zeros_like(prev_p) 
        dp[2:], _ = get_coeffs(beta, alpha) 
        new_p = prev_p + dp
        # recalculate A based on new a,b,c
        new_p[0] = get_A(x_data, Nsat, *new_p[2:])

        new_metric = metric_func(x_data, y_data, bins, sigma, func, *new_p)

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
        if j > 200:
            print('max iterations reached')
            return p, j
        if abs(delta_metric/len(initial_guess)) < acc:
            print(f'accuracy threshold reached in {j} iterations')
            return p, j

def fit_data_set(filename, initial_guess, acc=0.000001):
    xmin, xmax = 1e-4, 5
    radius, nhalo = readfile(filename)
    N_per_halo, edges, Nsat = bin_data(radius, nhalo, xmin, xmax)
    edges1 = edges.copy()
    bin_centers = (edges[:-1] + edges[1:]) * 0.5
    
    # create the input parameters.
    A = get_A(bin_centers, Nsat, *initial_guess)
    params_i = [A, Nsat, *initial_guess]
    sigma = np.sqrt(N_per_halo)

    # remove bins that are zero. We're throwing away data here, but i couldn't figure out how else to handle the zero bins.
    nonzero_mask = N_per_halo > 0
    bin_centers_clipped = bin_centers[nonzero_mask]
    N_per_halo_clipped = N_per_halo[nonzero_mask]
    sigma_clipped = sigma[nonzero_mask]
    edges_clipped = edges[np.append(nonzero_mask, True)]

    # the actual fitting
    p, j = levenberg_marquardt(bin_centers_clipped,
                             N_per_halo_clipped,
                             edges_clipped,
                             params_i, sigma_clipped, N_integrated, chi_squared, chi, 
                             w=20,
                             lamb = 1e-4,
                             acc=acc)
    edges = edges_to_bins(edges_clipped)

    # See what the minimum chi squared is
    chi2_min = chi_squared(bin_centers_clipped, N_per_halo_clipped, edges, sigma_clipped, N_integrated, *p)

    return Nsat, N_per_halo, p, chi2_min, bin_centers, edges1

# Plot of binned data with the best fit (question 1b.4 and 1c)

fig1b, ax = plt.subplots(3,2,figsize=(6.4,8.0))
datalist = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']
ylims = [(1e-5, 1), (1e-4, 1e1), (1e-4, 1e2), (1e-3, 1e2), (1e-1, 1e4)]
p_list_b = np.zeros((5, 5))
min_chi2_list_b = np.zeros(5)
for i in range(5):
    init_guess = [2, 0.8, 1.6]
    Nsat, N_per_halo, p, chi2_min, bin_centers, edges1 = fit_data_set(datalist[i], init_guess, acc=0.0001)
    edges = edges_to_bins(edges1)
    A_initial = get_A(bin_centers, Nsat, *init_guess)*2
    p_initial = [A_initial, Nsat, *init_guess]
    row=i//2
    col=i%2

    p_list_b[i] = p # save em for later
    min_chi2_list_b[i] = chi2_min
    
    
    ax[row,col].step(edges1[:-1], N_per_halo, where='post', label='binned data')
    ax[row,col].step(edges1[:-1], N_integrated(edges, *p), where='post', label='best-fit integrated N')
    ax[row,col].set(yscale='log', xscale='log', ylim = ylims[i], xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
    print(f'------------------------------\nFor file {datalist[i]}:')
    print('Using Levenberg-Marquardt to fit a chi squared')
    print(f'N_sat = ', Nsat)
    print(f'a = {p[2]}\nb = {p[3]}\nc = {p[4]}')
    print(f'A = {p[0]}')
    print(f'Minimum chi squared value = {chi2_min}')
    print(f'Minimum chi squared corrected for degrees of freedom = {chi2_min / (len(bin_centers) - len(p))}')

ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1b.png', dpi=600)

# 1c. Let's do an MCMC.
def xor_shift(seed, a1 = 21, a2 = 35, a3 = 4):
    '''Performs 64bit xor shift on the seed number.'''
    state = np.uint64(seed)
    a1 = np.uint64(a1)
    a2 = np.uint64(a2)
    a3 = np.uint64(a3)

    state = state ^ (state >> a1)
    state = state ^ (state << a2)
    state = state ^ (state >> a3)
    return state

def mult_w_carry(seed, a = 4294957665):
    '''Performs multiply with carry, base 2**32'''
    state = np.uint64(seed)
    a = np.uint64(a)

    state = a*(state & np.uint64(2**32 - 1)) + (state >> np.uint64(32))
    number = state & np.uint64(2**32 - 1)

    return number, state

def rng_int(seed, size):
    state = np.uint64(seed)
    result = np.zeros(size)
    for i in range(size):
        xor_state = xor_shift(state)
        if xor_state < 2**32:
            number, mult_state = mult_w_carry(xor_state)
        else:
            number, mult_state = mult_w_carry(xor_state & np.uint64(2**32 - 1))
        state = mult_state
        result[i] = number

    return result

def rng_float(seed, size, a, b):
    ints = rng_int(seed, size)
    floats = a + (b-a)*(ints / (2**32))
    return floats

def box_muller(x1, x2):
    r = np.sqrt(-2 * np.log(x1))
    theta = 2 * np.pi * x2
    z1 = r * np.cos(theta)
    z2 = r * np.sin(theta)
    return z1, z2

def multivariate_normal(mean, sigma):
    # instead of sigma, there should be the covariance matrix, but since i don't know anything
    # about the parameters, i'm just going to say they have variance of 1.
    # there is probably a nice way to do this, but i just don't know.
    n_walkers, n_params = mean.shape
    
    # get independent random samples
    z = np.zeros((n_walkers, n_params))
    u1_list = rng_float(i+413, n_walkers*n_params, 0, 1).reshape(n_walkers, n_params)
    u2_list = rng_float(i+13, n_walkers*n_params, 0, 1).reshape(n_walkers, n_params)
    for j in range(0, n_params, 2):
        u1 = u1_list[:, j]
        u2 = u2_list[:, j]
        z1, z2 = box_muller(u1, u2)
        z[:, j] = z1
        if j + 1 < n_params:
            z[:, j+1] = z2
    # go back to multivariate normal 
    samples = mean + sigma * z
    return samples

def q(p):
    sigma = np.ones_like(p) * 0.1
    new_p = multivariate_normal(p, sigma)
    # ensure Nsat stays the same
    new_p[:, 1] = p[:, 1]
    return new_p

def lnL(ydata, edges, func, p):
    lambda_i = np.array([np.array(func(edges, *params)) for params in p])
    lambda_i[lambda_i <= 0] = 1e-40 # make sure no invalid values are entered into the log
    ydata[ydata == 0] = 1e-40
    likelihood = ydata * np.log(lambda_i) - lambda_i
    lnL = np.sum(likelihood, axis = 1)
    return -lnL

def lnP(p):
    # log prior likelihood for a parameter p. let's do a flat prior, where
    # the values can't be negative.
    ln_P = np.zeros(p.shape[0])
    ln_P[np.any(p < 0, axis=1)] = -np.inf 
    return ln_P

def metropolis_hastings(ydata, edges, func, p_start, lnL, lnP, n_walkers, n_steps, filename):
    n_param = len(p_start)
    edges1 = edges_to_bins(edges)

    p = np.zeros((n_walkers, n_steps, n_param))
    p_init = np.linspace(p_start * 0.8, p_start * 1.2, n_walkers)
    
    p[:, 0, :] = p_init
    p[:, 0, 1] = np.ones_like(p[:,0,1])*p_start[1]

    y_list = rng_float(173937241, n_steps*n_walkers, 0, 1).reshape(n_walkers, n_steps)
    acceptance_list = np.zeros(n_steps)                                          
    for i in range(n_steps-1):
        p_prop = q(p[:, i, :])

        lnpi_i = lnL(ydata, edges1, func, p[:, i, :]) + lnP(p[:, i, :])  
        lnpi_prop = lnL(ydata, edges1, func, np.abs(p_prop)) + lnP(p_prop)

        log_L_ratio = lnpi_prop-lnpi_i
        y = y_list[:, i]
        
        accept = log_L_ratio >= 0.5*np.log(y)
        acceptance_list[i] = np.sum(accept)/len(accept)
        p[:, i + 1, :] = np.where(accept[:, None], p_prop, p[:, i, :]) 

        if i % 100 == 0:
            print(f'now at step {i}. {n_steps-i} remaining.')

    print('acceptance rate:', np.mean(acceptance_list))

    # discard the burn-in:
    p_full = p.copy()
    p = p[:, 100:, :]

    return p, p_full

def plot_mcmc(p, filename='my_mcmc_result.pdf'):
    n_param = 5
    param_names = ['A', 'Nsat', 'a', 'b', 'c']
    fig, ax1 = plt.subplots(n_param, 1,figsize=(6.4,8.0))
    for i in range(n_param):
        for j in range(p.shape[0]):
            ax1[i].plot(p[j, :, i], alpha=0.2, color='black')
            ax1[i].set(xlabel='step', ylabel=f'{param_names[i]}')
    plt.tight_layout()
    plt.savefig(filename, dpi=600)
    plt.close()

def fit_data_set_mcmc(filename, initial_guess):
    print('-----------------------')
    print('fitting data set', filename,'...')

    xmin, xmax = 1e-4, 5
    radius, nhalo = readfile(filename)
    N_per_halo, edges, Nsat = bin_data(radius, nhalo, xmin, xmax)
    edges1 = edges.copy()
    bin_centers = (edges[:-1] + edges[1:]) * 0.5
    
    # create the input parameters.
    A = get_A(bin_centers, Nsat, *initial_guess)
    params_i = np.array([A, Nsat, *initial_guess])

    # remove bins that are zero
    nonzero_mask = N_per_halo > 0
    N_per_halo_clipped = N_per_halo[nonzero_mask]
    edges_clipped = edges[np.append(nonzero_mask, True)]

    # the actual fitting
    p, p_full = metropolis_hastings(N_per_halo_clipped, edges_clipped, N_integrated, params_i,
                           lnL, lnP, n_walkers = 50, n_steps = 300, filename = filename)
    
    p_avg = np.mean(p, axis=(0, 1))
    edges = edges_to_bins(edges_clipped)

    # check the minimum log likelihood value
    lnL_min = lnL(N_per_halo_clipped, edges, N_integrated, np.array([p_avg]))
    
    return N_per_halo, p_avg, lnL_min, bin_centers, edges1, p_full

# Plot 1c (same code as above)
fig1b, ax = plt.subplots(3,2,figsize=(6.4,10.0))
datalist = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']
ylims = [(1e-5, 1), (1e-4, 1e1), (1e-4, 1e2), (1e-3, 1e2), (1e-1, 1e4)]

p_list_c = np.zeros((5, 5))
min_logL_list_c = np.zeros(5)
p_full_list = []
for i in range(5):
    init_guess = np.array([1.6, 0.9, 2.5])
    N_per_halo, p, lnL_min, bin_centers, edges1, p_full = fit_data_set_mcmc(datalist[i], init_guess)
    edges = edges_to_bins(edges1)
    A_initial = get_A(bin_centers, Nsat, *init_guess)
    p_initial = [A_initial, Nsat, *init_guess]
    row=i//2
    col=i%2
    
    # export some values
    p_list_c[i] = p
    min_logL_list_c[i] = lnL_min
    p_full_list.append(p_full)
    ax[row,col].step(edges1[:-1], N_per_halo, where='post', label='binned data')
    ax[row,col].step(edges1[:-1], N_integrated(edges, *p), where='post', label='best-fit integrated N')
    ax[row,col].set(yscale='log', xscale='log', ylim = ylims[i], xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
    print(f'------------------------------\nFor file {datalist[i]}:')
    print(f'Using MCMC to fit a Poisson likelihood')
    print(f'a = {p[2]}\nb = {p[3]}\nc = {p[4]}')
    print(f'A = {p[0]}')
    print(f'Minimum log likelihood value = {lnL_min}')

ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1c.png', dpi=600)

for i in range(5):
    plot_mcmc(p_full_list[i], filename=f'{datalist[i]}.pdf')
# 1D...... help

# we need continuous E, and discrete O. Undo the normalization.

def N_integrated_unnormalized(edges, A, Nsat, a, b, c, nhalo):
    result = np.zeros_like(edges[:,0])
    for i in range(len(result)):
        xmin, xmax = edges[i]
        result[i] = romberg_interval(N, xmin, xmax, 5, A, Nsat, a, b, c)[0] # * nhalo
    return result 

def bin_data_unnormalized(radius, nhalo, xmin, xmax):
    edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 51)
    binned_data, edges = np.histogram(radius, bins=edges)
    N_binned = binned_data * np.diff(edges)# don't normalize now

    # average Nsat
    Nsat = len(radius) / nhalo
    return N_binned, edges, Nsat

def G_test(O, E):
    return 2*np.sum(O * np.log(O) - np.log(E))

def chi2_dist(x, k):
    from scipy.special import gamma
    return x**(k/2 -1) * np.exp(-x/2) / (2**(k/2) * gamma(k/2))

def chi2_cdf(x, k):
    from scipy.special import gamma, gammainc
    return gammainc(k/2, x/2) / gamma(k/2)

for i in range(5):
    p_b = p_list_b[i]
    min_chi2_b = min_chi2_list_b[i]
    p_c = p_list_c[i]
    min_logL = min_logL_list_c[i]

    radius, nhalo = readfile(datalist[i])
    O, edges, Nsat = bin_data_unnormalized(radius, nhalo, 1e-4, 5)

    # remove bins that are zero
    nonzero_mask = O > 0
    O_clipped = O[nonzero_mask]
    edges_clipped = edges[np.append(nonzero_mask, True)]
    edges1 = edges_to_bins(edges_clipped)

    k = len(edges[:-1]) - 3 # three fitted parameters: a, b, c
    E_b = N_integrated_unnormalized(edges1, *p_b, nhalo)
    E_c = N_integrated_unnormalized(edges1, *p_c, nhalo)

    G_b = G_test(O_clipped, E_b)
    G_c = G_test(O_clipped, E_c)

    Q_b = 1 - chi2_cdf(G_b, k)
    Q_c = 1 - chi2_cdf(G_c, k)

    print('-------------------------')
    print(f'for file {datalist[i]}:')
    print(f'chi2: G = {G_b}, Q = {Q_b}')
    print(f'lnL: G = {G_c}, Q = {Q_c}')




