import numpy as np
import matplotlib.pyplot as plt
import copy 
def selectionsort(a):
    N = len(a)
    for i in range(0, N-1):
        imin = i
        for j in range(i+1, N):
            if a[j] < a[imin]:
                imin = j
        if imin != i:
            a[[imin, i]] = a[[i, imin]]
    return a

rand_arr = np.random.randint(0, 100, size = 10)

print('using selectionsort:')
print(rand_arr)
sorted_arr1 = selectionsort(rand_arr)
print(sorted_arr1)

# now implement mergesort.

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

print('using mergesort:')
rand_arr2 = np.random.randint(0, 100, size = 10)
print(rand_arr2)
sort_idx2 = mergesort(rand_arr2)
sorted_arr2 = rand_arr2[sort_idx2]
print(sorted_arr2)
print(sort_idx2)

def bisection(func, brac, *args, abs_tol = 1e-8, frac_tol = 1e-8, max_iter = 50):
    '''implements the bisection method of root finding.'''
    i = 0
    a = brac[0]
    b = brac[1]
    # check if a, b is a bracket
    assert func(a, *args)*func(b, *args) < 0

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

def secant(func, guesses, *args, abs_tol = 1e-8, frac_tol = 1e-8, max_iter = 50):
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

def false_position(func, brac, *args, abs_tol = 1e-8, frac_tol = 1e-8, max_iter = 50):
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

def x_linear(x):
    return 2*x - 4

def x_squared(x):
    return 2*x*x - 4 
brac = [-1, 8]
root, num_it = bisection(x_linear, brac)
print('using bisection:', root, num_it)

guesses = [500, 1000]
root2, num_it2 = secant(x_squared, guesses= guesses)
print('using secant:', root2, num_it2)

root3, num_it3 = false_position(x_squared, brac, max_iter = 1000)
print('using false position:', root3, num_it3)



