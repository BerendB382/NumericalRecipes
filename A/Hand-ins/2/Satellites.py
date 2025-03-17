import numpy as np
import matplotlib.pyplot as plt
import copy 

def n(x, A, Nsat, a, b, c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def N(x, A, Nsat, a, b, c):
    return n(x, A, Nsat, a, b, c)*4*np.pi*x**2

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

A=1. # to be computed
Nsat=100
a=2.4
b=0.25
c=1.6

# We want to integrate the function (with A = 1) to find what we need to set A to 
# to find <Nsat> = 100. A will be 100 divided by the value of the integral.

# Because n(x) is radially symmetric, we can use the fact that 
# dV = 4pi*r**2 dr to write down the integral.

x = np.linspace(0.00001, 5, 15)
integ_result = romberg(N, x, 10, A, 100, a, b, c)
new_A = Nsat / integ_result
print('A found by integral normalisation:', new_A)

check = romberg(N, x, 10, new_A, Nsat, a, b, c)
check = Nsat / check
print('Sanity check:', check)

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

def mergesort(a):
    N = len(a)
    idx = np.arange(0, N, 1)
    def mergesort_helper(idx):
        if len(idx) > 1:
            mid_idx = len(idx)//2
            left = copy.deepcopy(idx[:mid_idx])
            right = copy.deepcopy(idx[mid_idx:])

            # recursively sort the left and right arrays 
            mergesort_helper(left) 
            mergesort_helper(right)
            i, j, k = 0, 0, 0
            # go through the arrays and check which element of is the smallest between them, then add that one to a
            while i < len(left) and j < len(right):
                if a[left[i]] < a[right[j]]:
                    idx[k] = left[i]
                    i += 1
                else:
                    idx[k] = right[j]
                    j += 1
                k += 1
            # if left or right still hold elements after one of the lists is ran through, put them into a in order
            while i < len(left):
                idx[k] = left[i]
                i += 1
                k += 1
            while j < len(right):
                idx[k] = right[j]
                j += 1
                k += 1
    
    mergesort_helper(idx)
    return idx


#Plot of histogram in log-log space with line (question 1b)
xmin, xmax = 10**-4, 5
N_generate = 10000

# Make the analytical function
relative_radius = np.linspace(xmin, xmax, 10000) 
analytical_function = N(relative_radius, new_A, Nsat, a, b, c) / Nsat

# find the maximum value by sorting
sorting_idx2 = mergesort(analytical_function)
sorted_analytical = analytical_function[sorting_idx2]
max_val2 = sorted_analytical[-1]

# Do rejection sampling for r.
Pu3 = 10**rng_float(44, N_generate*10, np.log10(xmin), np.log10(xmax))
Pu4 = rng_float(33, N_generate*10, 0, max_val2)
P_x = N(Pu3, new_A, Nsat, a, b, c) / Nsat 

rej_result = []
for i in range(N_generate*10):
    if Pu4[i] < P_x[i]:
        rej_result.append(Pu3[i])
        if len(rej_result) == 10000:
            break
rej_result = np.array(rej_result) # take the first 10000 elements

#21 edges of 20 bins in log-space
edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 21)
hist = np.histogram(rej_result, bins=edges)[0]
hist_scaled = Nsat * hist / (len(Pu3)) # normalize by the amount of samples originally used. 

hist_px = np.histogram(P_x, bins = edges)[0]
hist_px_scaled = Nsat * hist / (len(P_x))

fig1b, ax = plt.subplots()
# ax.stairs(hist_px_scaled, edges = edges, fill=False, label = 'calculated Px')
ax.stairs(hist_scaled, edges=edges, fill=False, label='Satellite galaxies') 
plt.plot(relative_radius, analytical_function, 'r-', label='Analytical solution') 
ax.set(xlim=(xmin, xmax), ylim=(10**(-3), 10), yscale='log', xscale='log',
       xlabel='Relative radius', ylabel='Number of galaxies')
ax.legend()
plt.savefig('my_solution_1b.png', dpi=600)
plt.close()

# We'll achieve this by generating a random sequence using our random number generator.
# We will then sort the array and thus generate a shuffle index which we then apply to
# the samples from above.
rnd_seq = rng_float(69, 10000, 0, 1)
rnd_shuffle_idx = mergesort(rnd_seq)

# Select 100 random satellite galaxies by applying the index to the samples
# and then selecting the first 100
chosen = rej_result[rnd_shuffle_idx][:100]
sorting_idx = mergesort(chosen)
chosen_sorted = chosen[sorting_idx]

fig1c, ax = plt.subplots()
ax.plot(chosen_sorted, np.arange(100))
ax.set(xscale='log', xlabel='Relative radius', 
       ylabel='Cumulative number of galaxies',
       xlim=(xmin, xmax), ylim=(0, 100))
plt.savefig('my_solution_1c.png', dpi=600)

# Now we will calculate dn(x)/dx at x = 1. 

def n_deriv(x, A, Nsat, a, b, c):
    return A*Nsat/b * ((a-3)*(x/b)**(a-4) - c*(x/b)**(a+c-4))*np.exp(-(x/b)**c)

def central_diff(x, func, h, *args):
    deriv = (func(x + h, *args) - func(x - h, *args)) / (2*h)
    return deriv

def ridder(x, func, h, d, *args, order = 5, tol = 1e-8):
    r = np.zeros(order)
    oneDeeth = 1 / d
    r[0] = central_diff(x, func, h, *args)
    for o in range(1, order):
        h = h * oneDeeth
        r[o] = central_diff(x, func, h, *args)
    # now iterate to combine the previous 
    
    for j in range(1, order):
        for i in range(order-1-j):
            r[i] = (d**(2*j) * r[i+1] - r[i]) / (d**(2*j) - 1)
    return r[0]

#Cumulative plot of the chosen galaxies (1c)
x = 1 
h = 0.2

deriv_x = ridder(x, n, h, 2, new_A, Nsat, a, b, c, order = 10)
deriv_x_a = n_deriv(x, new_A, Nsat, a, b, c)
print('derivative of n with respect to x')
print("with Ridder's method:", deriv_x)
print('Analytically:', deriv_x_a)