import numpy as np
import matplotlib.pyplot as plt

def FFT1D(a):
    N = np.int64(len(a))
    if (N & (N-1) != 0) and N == 0:
        raise Exception('N is not a power of two. Please input an array that is a power of two long.')
    def FFT1D_helper(a, N):
        a = np.array(a, dtype=np.complex64)
        Nth = 1/N
        N = int(N)
        if N > 2:
            even_a = a[::2]
            odd_a = a[1::2]

            even_a = FFT1D_helper(even_a, N*0.5)
            odd_a = FFT1D_helper(odd_a, N*0.5)

            a = np.concatenate((even_a, odd_a))
        for k in range(int(N*0.5)):
            coolexp = np.exp(2j*np.pi*k*Nth)
            t = a[k]
            a[k] = t + coolexp*a[int(k+N*0.5)]
            a[int(k+N*0.5)] = t - coolexp*a[int(k+N*0.5)]
        return a
    a = FFT1D_helper(a, N)
    return a

def inv_FFT1D(a):
    N = np.int64(len(a))
    if (N & (N-1) != 0) and N == 0:
        raise Exception('N is not a power of two. Please input an array that is a power of two long.')
    a = a
    def FFT1D_helper(a, N):
        a = np.array(a, dtype=np.complex64)
        Nth = 1/N
        N = int(N)
        if N > 2:
            even_a = a[::2]
            odd_a = a[1::2]

            even_a = FFT1D_helper(even_a, N*0.5)
            odd_a = FFT1D_helper(odd_a, N*0.5)

            a = np.concatenate((even_a, odd_a))
        for k in range(int(N*0.5)):
            coolexp = np.exp(2j*np.pi*k*Nth)
            t = a[k]
            a[k] = t - coolexp*a[int(k+N*0.5)]
            a[int(k+N*0.5)] = t + coolexp*a[int(k+N*0.5)]
        return a
    a = FFT1D_helper(a, N)
    return a

# notities voor toekomstige berend:
# je IFFT gaat nog mis met een soort schaalfactor en also hoe schijt het is. Dubbelcheck of je FFT klopt
# en ga dan verder met je IFFT! TODO: TODO TODO TODO 
def f(x):
    return np.complex64((2*x + np.sin(2*np.pi*x*0.2) + 3*np.cos(2*np.pi*x*0.5))*np.sin(2*x))

def f(x):
    return np.sin(2*np.pi*5*x)

x = np.linspace(0, 1, 1024)
arr = f(x)
# We have to shift the output of our function, because we shifted it in calculation to make it faster.
# put the right half before the first half.
def putback(F):
    return np.concatenate((F[int(len(F)/2):], F[:int(len(F)/2)]))
F = FFT1D(arr)
fouriers_shifted = putback(F)

k = np.arange(-len(F)/2, len(F)/2, 1)

plt.figure(figsize=(12, 8))
# plt.plot(k, abs(fouriers), label='FFT Coeff amplitude')
plt.stem(k, abs(fouriers_shifted), label='FFT Coeff amplitude')
plt.xlabel('Fourier space variable k')
plt.ylabel('FFT coefficients')
plt.title('FFT results')
plt.legend()
plt.show()
plt.close()

F_back = inv_FFT1D(F)

plt.plot(x, arr, label='original signal')
plt.plot(x, F_back, label='reconstructed function')
plt.xlabel('x')
plt.ylabel('y')
plt.title('reconstruction')
plt.legend()
plt.show()
