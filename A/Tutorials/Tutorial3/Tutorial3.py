import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.linalg import lu

def bisection2(x, x_values):
    low, high = 0, len(x_values) - 1
    while high - low > 1:
        mid = (low + high) // 2
        if x_values[mid] > x:
            high = mid
        else:
            low = mid
    return high, low

def linear_interp(x_new, x_values, y_values):
    y_new = []
    delta_x = x_values[5] - x_values[4]
    for x in x_new:
            i, j = bisection2(x, x_values)
            a = (int(y_values[i]) - int(y_values[j]))
            y_interp = y_values[j] + (x - x_values[j]) * a 
            y_new.append(y_interp)
    return np.array(y_new)

# image = imread('/Users/bjhnieuwhof/Google Drive/Universiteit Leiden/Master Astronomy/Numerical Recipes/A/Tutorials/Tutorial1+2/M42_128.jpg')
# firstline = image[0]

# old_x = np.linspace(0, len(firstline), len(firstline))
# new_x = np.linspace(0, len(firstline), 201)
# new_y = linear_interp(new_x, old_x, firstline)

# plt.plot(old_x, firstline, label = 'old line')
# plt.plot(new_x, new_y, label = 'interpolated line')
# plt.xlabel('x')
# plt.xlabel('y')
# plt.title('linear interpolation result')
# plt.legend()
# plt.show()

def crout(A):
    '''Overwrites a matrix by its LU composition.'''
    assert A.shape == A.T.shape
    alpha = np.identity(len(A))
    for j in range(len(A)):
        for i in range(len(A)):
            if i <= j:
                if i >= 0:
                    sumvar_a = np.arange(0, i, 1)
                    sum_a = np.sum(alpha[i, sumvar_a]*A[sumvar_a, j])
                    A[i, j] -= sum_a
                
            else:
                if j > 0:
                    sumvar_b = np.arange(0, j, 1)
                else:
                    sumvar_b = np.array([int(0)]) 
                sum_b = np.sum(alpha[i, sumvar_b]*A[sumvar_b, j])
                A[i, j] = 1 / A[j, j] * (A[i, j] - sum_b)
                if i != j:
                    alpha[i, j] = A[i, j]
        
    return A

A = np.array([[1, -4, 2],
              [-2, 1, 3],
              [2, 6, 8]])

print(A)

print(crout(A))

# P, L, U = lu(A)

# print('P scipy:\n', P)
# print('L scipy:\n', L)
# print('U scipy:\n', U)
# print('scipy check:\n', P @ L @ U)


        




