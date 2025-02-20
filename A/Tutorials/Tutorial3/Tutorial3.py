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

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.image as mpimg

img = mpimg.imread('A/Tutorials/Tutorial3/petiti_luigi.jpg')
grey_img = np.dot(img[..., :3], [0.640, 0.595, 0.155]) # Rec601 hack to turn greyscale

U, S, Vh = np.linalg.svd(grey_img, full_matrices=False)
print(U.shape, np.diag(S).shape, Vh.shape)
remake_img = U @ np.diag(S) @ Vh

# plt.imshow(grey_img, cmap='gray')
# plt.show()
# plt.imshow(remake_img, cmap='gray')
# plt.show()

size_lim = 1000

U, S, Vh = U[:size_lim, :size_lim], S[:size_lim], Vh[:size_lim, :size_lim]
print(U.shape, np.diag(S).shape, Vh.shape)
compressed_img = U @ np.diag(S) @ Vh

# plt.imshow(compressed_img, cmap='gray')
# plt.show()

# EXERCISE 2
class Matrix():
    def __init__(self, n_rows, n_cols):
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.arr = np.zeros((n_rows, n_cols))

    def row_swap(self, i, j):
        if i != j:
            row_i = np.copy(self.arr[i, :])
            self.arr[i, :] = self.arr[j, :]
            self.arr[j, :] = row_i 
    
    def LU_comp(self):
        return crout(self.arr)
    
    def solve_sys_GJ(self):
        for col_idx in range(self.n_cols-1):
            pivot_idx = np.argmax(self.arr[:, col_idx])
            if math.isclose(self.arr[pivot_idx, col_idx], 0):
                break
            self.row_swap(pivot_idx, col_idx)
            self.arr[col_idx, :] = self.arr[col_idx, :]/self.arr[col_idx, col_idx]

            # use the pivot row to reduce all other rows with non-zero elements in column i.
            for row_idx in np.arange(col_idx+1, self.n_rows):
                if row_idx >= self.n_rows or col_idx >= self.n_cols:
                    break
                if not math.isclose(self.arr[row_idx, col_idx], 0):
                    self.arr[row_idx, :] -= self.arr[col_idx, :]*self.arr[row_idx, col_idx]
        return self.arr
                





A = Matrix(5, 5)
A.arr = np.array([[3, 8, 1, -12, -4, 2],
                  [1, 0, 0, -1, 0, 0],
                  [4, 4, 3, -40, -3, 1],
                  [0, 2, 1, -3, -2, 0],
                  [0, 1, 0, -12, 0, 0]]).astype(np.float32)

gj_a = A.solve_sys_GJ()
print(np.array(gj_a))




        




