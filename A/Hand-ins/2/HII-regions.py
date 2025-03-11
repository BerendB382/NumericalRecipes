import numpy as np
import matplotlib.pyplot as plt
import copy 

k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s

# here no need for nH nor ne as they cancel out
def equilibrium1(T,Z,Tc,psi):
    return psi*Tc*k - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T*k


def equilibrium2(T,Z,Tc,psi, nH, A, xi):
    return (psi*Tc - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T - .54 * ( T/1e4 )**.37 * T)*k*nH*aB + A*xi + 8.9e-26 * (T/1e4)

def bisection(func, brac, *args, abs_tol, frac_tol, max_iter):
    i = 0
    a = brac[0]
    b = brac[1]
    while i < max_iter:
        # Divide the bracket in 2 
        c = (a+b)/2 
        if func(a, *args)*func(c, *args) < 0:
            b = c
        else:
            a = c
        i += 1
        if a-b < abs_tol or (a-b)/a < frac_tol:
            return a
    
    return a, i

