import numpy as np
import matplotlib.pyplot as plt

def trap_integrate(N, func, a, b): 
    h = (b - a)/N
    x = np.linspace(a, b, N)
    y = func(x)
    result = h * ((y[0]+y[-1])*0.5 + np.sum(y[1:-1]))
    return result

def x_squared(x):
    return x**2

def sin(x):
    return np.sin(x)

square_integrated = trap_integrate(2000, x_squared, 1, 5)

def simpson(N, func, a, b):
    s0 = trap_integrate(N, func, a, b)
    s1 = trap_integrate(2*N, func, a, b)

    s4 = 1/3*(4*s1 - s0)
    return s4

simp_result = simpson(2000, x_squared, 1, 5)
print(square_integrated)
print(simp_result)

def romberg(func, a, b, order):
    r = np.zeros(order)
    
    # calculate r0
    h = (b - a)/N
    r[0] = 0.5 * h * (func(a) - func(b))
    Np = 1
    for i in np.arange(1, order-1):
        r[i] = 0
        delta = h
        h = 0.5 * h
        x = a + h
        for n in range(Np):
           r[i] += func(x)
           x += delta
        r[i] = 0.5*(r[i-1] + delta*r[i])
        Np = 2*Np
    Np = 1
    for i in np.arange(1, order-1):
        Np = 4*Np 
        for j in np.arange(0, order-i):
            r[j] = (Np * r[j+1] - r[j]) / (Np - 1)

    return r

romberg_result = romberg(x_squared, 1, 5, 15)
x = np.linspace(1, 5, 20)
trap_result = np.trapezoid(x_squared(x), x)

print(romberg_result)
print(trap_result)

# je moet deze functie dus per interval gaan callen.