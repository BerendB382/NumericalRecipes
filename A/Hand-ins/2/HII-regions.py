import numpy as np
import matplotlib.pyplot as plt
import copy 
import timeit

def bisection(func, brac, *args, abs_tol = 1e-8, frac_tol = 1e-8, max_iter = 50):
    '''implements the bisection method of root finding.'''
    i = 0
    a = brac[0]
    b = brac[1]
    # check if a, b is a bracket
    # assert func(a, *args)*func(b, *args) < 0

    while i < max_iter:
        # Divide the bracket in 2 
        c = (a+b)/2 
        if func(a, *args)*func(c, *args) < 0:
            b = c
        else:
            a = c
        i += 1
        if abs(a-b) < abs_tol or abs((a-b)/a) < frac_tol:
            return a, i
    return a, i

def secant(func, guesses, *args, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 50):
    '''implements the secant method for root finding.'''
    i = 2
    # r will hold all guesses
    r = np.zeros(max_iter)
    r[0] = guesses[0]
    r[1] = guesses[1]
    while i < max_iter:
        r[i] = r[i-1] - ((r[i-1] - r[i-2])/(func(r[i-1], *args) - func(r[i-2], *args))) * func(r[i-1], *args)
        if abs(r[i] - r[i-1]) < abs_tol or abs((r[i] - r[i-1])/r[i]) < frac_tol:
            return r[i], i
        i += 1
    
    print('maximum secant interations reached')
    return r[i-1], i

def false_position(func, brac, *args, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 50):
    '''implements the false position method for root finding.'''
    i = 2
    r = np.zeros(max_iter)
    r[0] = brac[0]
    r[1] = brac[1]
    # check if a, b is a bracket
    assert func(r[0], *args)*func(r[1], *args) < 0

    while i < max_iter:
        r[i] = r[i-1] - ((r[i-1] - r[i-2])/(func(r[i-1], *args) - func(r[i-2], *args))) * func(r[i-1], *args)
        if func(r[i], *args)*func(r[i-2], *args) < 0:
            r[i-1] = r[i-2]
        if abs(r[i] - r[i-1]) < abs_tol or abs((r[i] - r[i-1])/r[i]) < frac_tol:
            return r[i], i
        i += 1
    
    return r[i-1], i

## 2a. 

k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
Z = 0.015 # unitless
psi = 0.929 # unitless
Tc = 1e4 # K

# here no need for nH nor ne as they cancel out
def equilibrium1(T,Z,Tc,psi):
    return psi*Tc*k - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T*k 

# test and visualize first
temps = np.linspace(1, 1e7, 1000)
eqs = equilibrium1(temps, Z, Tc, psi)
# plt.plot(temps, eqs)
# plt.show()

brac1 = [1, 1e7] # K 
print(brac1)
midpoint1 = 0.5e7 # mid point of the bracket interval
print('------------2a--------------')

root1, num_it1 = bisection(equilibrium1, brac1, Z, Tc, psi, abs_tol = 0.01, frac_tol = 1e-8, max_iter = 50)
print('root and number of iterations')
print('using bisection:', root1, num_it1)
print('function evaluated at root:', equilibrium1(root1, Z, Tc, psi))

root2, num_it2 = secant(equilibrium1, brac1, Z, Tc, psi, abs_tol = 0.01, frac_tol = 1e-8, max_iter = 50)
print('using secant method:', root2, num_it2)
print('function evaluated at root:', equilibrium1(root2, Z, Tc, psi))

root3, num_it3 = false_position(equilibrium1, brac1, Z, Tc, psi, abs_tol = 0.01, frac_tol = 1e-8, max_iter = 50)
print('using false position method:,', root3, num_it3)
print('function evaluated at root:', equilibrium1(root3, Z, Tc, psi))

def equilibrium2(T,Z,Tc,psi, nH, A, xi):
    return (psi*Tc - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T - .54 * ( T/1e4 )**.37 * T)*k*nH*aB + A*xi + 8.9e-26 * (T/1e4)


print('---------------2b----------------')
# 2b
A = 5e-10 # erg
xi = 1e-15 # s**-1
brac2 = [1, 1e15] # K
k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
Z = 0.015 # unitless
psi = 0.929 # unitless
Tc = 1e4 # K

# test and visualize first
temps = np.linspace(1, 1e15, 1000)
ne = 1e5
eqs = equilibrium2(temps, Z, Tc, psi, ne, A, xi)
# plt.plot(temps, eqs)
# plt.show()

for ne in [1e-4, 1, 1e4]:
    print('\nfor ne = ', ne)
    root1, num_it1 = bisection(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100)
    print('root and number of iterations')
    print('using bisection:', root1, num_it1)
    print('function evaluated at root:', equilibrium2(root1, Z, Tc, psi, ne, A, xi))

    root2, num_it2 = secant(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100)
    print('using secant method:', root2, num_it2)
    print('function evaluated at root:', equilibrium2(root2, Z, Tc, psi, ne, A, xi))

    root3, num_it3 = false_position(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100000)
    print('using false position method:,', root3, num_it3)
    print('function evaluated at root:', equilibrium2(root3, Z, Tc, psi, ne, A, xi))

print('-------------------------------')
def smallest_ne():
    ne = 1e-4
    root1, num_it1 = bisection(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100)
    return root1, num_it1

def middle_ne():
    ne = 1
    root2, num_it2 = secant(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100)
    return root2, num_it2
def big_ne():
    ne = 1e4
    root2, num_it2 = secant(equilibrium2, brac2, Z, Tc, psi, ne, A, xi, abs_tol = 1e-10, frac_tol = 1e-10, max_iter = 100)
    return root2, num_it2

small_ne_time = timeit.timeit(lambda : smallest_ne(), number = 1)
mid_ne_time = timeit.timeit(lambda : middle_ne(), number = 1)
big_ne_time = timeit.timeit(lambda : big_ne(), number = 1)

print('Bisection for smallest n took:', small_ne_time)
print('Secant for middle ne took:', mid_ne_time)
print('Secant for biggest ne took:', big_ne_time)

temps = np.linspace(1, 1e15, 1000)
ne = 1e-4
eqs = equilibrium2(temps, Z, Tc, psi, ne, A, xi)
plt.grid()
plt.vlines(smallest_ne()[0], -4e-15, 1e-15, color='black')
plt.plot(temps, eqs)
plt.show()


