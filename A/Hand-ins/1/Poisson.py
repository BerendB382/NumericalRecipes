'''This file implements a function that returns the poisson distribution for an integer k.'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson



def log_factorial(k):
    '''Computes the factorial of a number in log space.'''
    if k == 0:
        result = 0
    else:
        k_values = np.arange(k, 0, -1)
        result = np.sum(np.log(k_values))
    return np.float32(result)

def log_poisson(lamb, k):
    '''Computes the log of the poisson distribution for a positive mean lambda.'''
    assert lamb > 0
    k = np.int32(k)

    result = k * np.log(lamb) - lamb - log_factorial(k)
    return np.float32(result)

def poisson(lamb, k):
    print(log_poisson(lamb, k))
    return np.float32(np.exp(log_poisson(lamb, k)))

lamb_k_values = np.array([np.array([1, 5, 3, 2.6, 100, 101], dtype=np.float32),
                          np.array([0, 10, 21, 40, 5, 200], dtype=np.int32)])

np.set_printoptions(precision=8)

for i in range(len(lamb_k_values[0])):
    lamb = lamb_k_values[0, i]
    k = lamb_k_values[1, i]
    poisson_val = poisson(lamb, k)
    print(f'for lamb = {lamb} and k = {k}, P_lambda(k) = {poisson_val}')


        




    
