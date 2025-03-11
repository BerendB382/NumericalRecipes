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

