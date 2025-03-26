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

def chi_squared(x_data, y_data, sigma, model, *args):
    return np.sum((y_data - model(x_data, *args))**2 / sigma**2)

def f(x, a, b, c):
    return a/x + b*x + c

def f_deriv_a(x):
    return 1/x

def f_deriv_b(x):
    return x

def f_deriv_c(x):
    return 1

def chi(x_data, y_data, sigma, model, deriv, *args):
    return np.sum(deriv(x_data)*(y_data - model(x_data, *args)) / sigma**2)

def calculate_alpha(x_data, p, sigma, lamb):
    n_param = len(p)
    alpha = np.zeros((n_param, n_param))
    funclist = [f_deriv_a, f_deriv_b, f_deriv_c]
    for k in range(n_param):
        for j in range(n_param):
            alpha[k, j] = np.sum((1/sigma**2)*(funclist[k](x_data)*funclist[j](x_data)))
    for i in range(n_param):
        alpha[i, i] = alpha[i, i] * (1+lamb)
    # TODO: somehow, your final element isn't correct. why?!?!
    return alpha

def calculate_beta(x_data, y_data, p, sigma, func):
    n_param = len(p)
    beta = np.zeros(n_param)
    xderivlist = [f_deriv_a, f_deriv_b, f_deriv_c]
    for k in range(n_param):
        beta[k] = chi(x_data, y_data, sigma, func, xderivlist[k], *p) # * -1
    
    return beta

def levenberg_marquardt(x_data, y_data, initial_guess, sigma, func, lamb=1e-3, w = 10):
    chi2 = chi_squared(x_data, y_data, sigma, func, *initial_guess)
    lamb = lamb
    j = 0
    p = initial_guess
    while True:
        alpha = calculate_alpha(x_data, p, sigma, lamb)
        beta = calculate_beta(x_data, y_data, p, sigma, func)
        prev_p = np.copy(p)
        dp, _ = get_coeffs(beta, alpha)
        new_p = prev_p + dp

        new_chi2 = chi_squared(x_data, y_data, sigma, func, *new_p)

        delta_chi = abs(new_chi2 - chi2)

        if new_chi2 >= chi2:
            lamb *= w 
            chi2 = chi2
            p = p
        else: 
            chi2 = new_chi2
            p = new_p
            lamb /= w
        j += 1
        if delta_chi/len(p) < 0.0001:
            return p, j
    
# a = 2
# b = 1
# c = 2

# x = np.linspace(0.5, 4, 20)

# xs = np.array([x for _ in range(1000)])

# data = f(xs, a, b, c)
# np.random.seed(42)
# gauss_noise = np.random.normal(0, 2, (1000, 20))

# noisy_data = data + gauss_noise

# plt.plot(xs[0], data[0], label='model data')
# plt.scatter(xs, noisy_data, color = 'red', label='noisy data')
# plt.title('')
# plt.xlabel('x')
# plt.ylabel('y')
# # plt.show()
# plt.close()


# plt.title('Many fits')
# plt.xlabel('x')
# plt.ylabel('y')

# initial_guess = [3, 2, 1]
# linny = np.linspace(0.5, 4, 50)
# params_arr = np.zeros((len(xs), len(initial_guess)))
# for i in range(len(xs)):
#     params, _ = levenberg_marquardt(xs[i], noisy_data[i], initial_guess, 2, f)
#     new_y = f(linny, *params)
#     plt.plot(linny, new_y, color='black', alpha=0.05)
#     params_arr[i] = params

# avg_params = np.mean(params_arr, axis = 0)
# avg_y = f(linny, *avg_params)
# plt.plot(xs[0], data[0], label = 'model data')
# plt.plot(linny, avg_y, linestyle = '--', color='red', label='average fit')
# plt.legend()
# plt.savefig('many_fits.pdf', dpi = 600)
    
smfdata = np.fromfile('A/Tutorials/Tutorial8/smfdata.txt', sep = '\n')
print(smfdata)

binned, edges = np.histogram(smfdata, 20, range=(8.5, 12.), density=True)

plt.stairs(binned, edges, fill = True, label = 'Stellar Mass Function')
plt.title('Stellar Mass Function')
plt.xlabel('log10(mass)')
plt.ylabel('counts')
plt.show()

        
