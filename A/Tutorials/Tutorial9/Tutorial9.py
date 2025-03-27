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

def z(x_data, y_data, sigma, model, *args):
    return (y_data - model(x_data, *args)) / sigma

def chi_squared(x_data, y_data, sigma, model, *args):
    return np.sum(z(x_data, y_data, sigma, model, *args)**2)

def lorentzian(x_data, y_data, sigma, model, *args):
    return np.sum(np.log(1 + 0.5 * z(x_data, y_data, sigma, model, *args)**2))

def lor(x_data, y_data, sigma, model, deriv, *args):
    return np.sum(deriv(x_data)/ sigma * z(x_data, y_data, sigma, model, *args) / (1 + 0.5*z(x_data, y_data, sigma, model, *args)**2))

def chi(x_data, y_data, sigma, model, deriv, *args):
    return np.sum(2*deriv(x_data)*(y_data - model(x_data, *args)) / sigma**2)

def f(x, a, b):
    return a*x + b

def f_deriv_a(x):
    return x

def f_deriv_b(x):
    return 1

def calculate_alpha(x_data, p, sigma, lamb):
    n_param = len(p)
    alpha = np.zeros((n_param, n_param))
    funclist = [f_deriv_a, f_deriv_b]
    for k in range(n_param):
        for j in range(n_param):
            alpha[k, j] = np.sum((funclist[k](x_data)*funclist[j](x_data)))
    for i in range(n_param):
        alpha[i, i] *= (1 + lamb)
    return alpha/sigma**2

def alpha_daan(x, sigma, lamb):
    func_list = [f_deriv_a, f_deriv_b]
    N = len(func_list)
    a = np.zeros([N,N])
    for i in range(len(x)):
        for j in range(N):
            for k in range(N):
                a[j,k]+=func_list[j](x[i])*func_list[k](x[i])
    a = a / sigma**2
    for i in range(len(a)):
            a[i, i] *= (1+lamb)
    return a

def calculate_beta(x_data, y_data, p, sigma, func, weight):
    n_param = len(p)
    beta = np.zeros(n_param)
    xderivlist = [f_deriv_a, f_deriv_b]
    for k in range(n_param):
        beta[k] = 0.5*weight(x_data, y_data, sigma, func, xderivlist[k], *p) # * -1
    return beta

def levenberg_marquardt(x_data, y_data, initial_guess, sigma, func, metric_func, weight_func, lamb=1e-3, w = 10, acc = 0.001):
    metric = metric_func(x_data, y_data, sigma, func, *initial_guess)
    lamb = lamb
    j = 0
    p = initial_guess
    while True:
        alpha = alpha_daan(x_data, sigma, lamb)
        print(alpha)
        alpha = calculate_alpha(x_data, p, sigma, lamb)
        print(alpha)
        beta = calculate_beta(x_data, y_data, p, sigma, func, weight_func)
        prev_p = np.copy(p)
        dp, _ = get_coeffs(beta, alpha)
        new_p = prev_p + dp

        new_metric = metric_func(x_data, y_data, sigma, func, *new_p)

        delta_metric = abs(new_metric - metric)
        print(delta_metric)
        if new_metric >= metric:
            lamb *= w 
            metric = metric
            p = p
        else: 
            metric = new_metric
            p = new_p
            lamb /= w
        j += 1
        if delta_metric/len(p) < acc:
            return p, j
    

datasets = np.zeros((5, 100, 2))
for i in range(1, 6):
    data = np.array(np.loadtxt('A/Tutorials/Tutorial9/outliers_dataset{0}.txt'.format(i)))
    datasets[i-1] = data

for i, data in enumerate(datasets):
    
    x = data[:, 0]
    y = data[:, 1]
    plt.scatter(x, y, alpha = 0.3, label = 'dataset {0}'.format(i))

plt.legend()
# plt.show()
plt.close()

data = datasets[0]
x = data[:, 0]
y = data[:, 1]
initial_guess = [14, 600]

plt.scatter(x, y, label='data')
plt.plot(x, f(x, *initial_guess), label='initial_guess', color='blue')
plt.xlabel('x')
plt.ylabel('y')

p_chi, num_it_chi = levenberg_marquardt(x, y, initial_guess, 1, f, chi_squared, chi, lamb = 0.001, acc = 110)
print('num iterations chi:', num_it_chi)
p_lor, num_it_lor = levenberg_marquardt(x, y, p_chi, 1, f, lorentzian, lor, lamb = 0.001, acc = 11)
print('num iterations lor:', num_it_lor)

plt.plot(x, f(x, *p_lor), label='fit w/ lorentzian', color='red')
plt.plot(x, f(x, *p_chi), label='fit w/ chi_squared', color = 'orange')

plt.legend()
plt.show()
plt.close()
