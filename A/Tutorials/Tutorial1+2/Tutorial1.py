import numpy as np
import matplotlib.pyplot as plt
import timeit

# %% 
def factorial(k):
    result = 1
    for i in range(k):
        result = result * k 
        k -= 1
    return result

def sinc_ps_vec(x, n_terms):
    terms = np.arange(0, n_terms, 1)
    terms = ((-1)**terms * x**(2*terms)) / factorial(2*terms+1)
    return np.sum(terms)

def sinc_ps_loop(x, n_terms):
    result = .0
    for n in range(n_terms):
        result += ((-1)**n * x**(2*n) ) / factorial(2*n+1)
    return result

def sinc_lib(x):
    return np.sin(x)/x

x = 7

# print('sinc_ps_vec:', sinc_ps_vec(x, 5))
print('sinc_lib:', sinc_lib(x))
print('sinc_ps_loop', sinc_ps_loop(x, 15))

n_terms_array = np.arange(0, 15, 1).astype(int)
print(n_terms_array)
error = np.zeros_like(n_terms_array)
for n in n_terms_array:
    error[n] = sinc_lib(x) - sinc_ps_loop(x, n)

plt.plot(error)
plt.xlabel('Number of terms in power series')
plt.ylabel('error')
plt.show()

# %%
def swr(M):
    return 2 * 6.67e-11 * M * 1.1111111e-17 


# %% 
from matplotlib.image import imread
import numpy as np
image = imread('/Users/bjhnieuwhof/Google Drive/Universiteit Leiden/Master Astronomy/Numerical Recipes/A/Werkcollege/M42_128.jpg')
first_line = image[0]

def lin_interp(j_low, j_high, x):
    '''interpolates linearly between j_low and j_high
    and provides the value at x'''
    a = (j_high[1] - j_low[1])/(j_high[0] - j_low[0])
    b = j_low[1] - a * j_low[0]
    return a*x + b

def bisection2(index, line):
    return [line[int(np.floor(index)), line[int(np.floor(np.ceil(index)))+1]]]
            

class interpolator():
    def __init__(self, data):
        self.data = data
    def interpolate_line(line, new_size):
        inter = np.zeros(new_size)
        indices = np.linspace(0, len(line)-1, new_size)
        for i in range(len(indices)):
            all_points = bisection2(i, line)
            a = all_points[1] - all_points[0]
            inter[i] = a*x+all_points[0]


        return inter

def bisection(data):
   return 0

# %%
interp = interpolator(image)
new_line = interpolator.interpolate_line(first_line, 201)

plt.plot(first_line)
plt.plot(new_line)
plt.title('linear interpolation result vs original line')
plt.show()

# %%
