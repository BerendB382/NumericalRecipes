import numpy as np
import matplotlib.pyplot as plt
import copy 

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
    
def brent(func, brac, *args, acc = 1e-5):
    return None

x = np.linspace(-10, 10, 1000)

# plt.show()
brac = [10, 20]
new_brac = find_bracket(f, brac)
print(new_brac)
minimum1 = golden_section(f, new_brac)
print(minimum1)
plt.plot(x, f(x))
plt.xlabel('x')
plt.ylabel('y')
plt.vlines(new_brac, 0, 5000)
plt.scatter(minimum1, f(minimum1), color='red')
plt.show()

