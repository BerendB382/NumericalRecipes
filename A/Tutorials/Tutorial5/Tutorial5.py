import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**2 * np.sin(x)
def f_deriv_a(x):
    return x**2 * np.cos(x) + 2*x * np.sin(x)

x = np.linspace(0, 2*np.pi, 1000)

y_func = f(x)
y_deriv_a = f_deriv_a(x)

def central_diff(x, func, h):
    deriv = (func(x + h) - func(x - h)) / (2*h)
    return deriv

y_deriv_cd = np.zeros((3, len(x)))
h_arr = np.array([0.1, 0.01, 0.001])
for h_idx, h in enumerate(h_arr):
    for i in range(len(x)):
        y_deriv_cd[h_idx, i] = central_diff(x[i], f, h) 
    # plt.plot(x, y_deriv_cd[h_idx], label=f'deriv w/ h = {h}') # normals
    plt.plot(x, y_deriv_a-y_deriv_cd[h_idx], label = f'deriv w/ h = {h}') # residuals


def ridder(x, func, h, d, order = 5):
    r = np.zeros(order)
    oneDeeth = 1 / d
    r[0] = central_diff(x, func, h)
    for o in range(1, order):
        h = h * oneDeeth
        r[o] = central_diff(x, func, h)
    # now iterate to combine the previous 
    for j in range(1, order):
        for i in range(order-1-j):
            r[i] = (d**(2*j) * r[i+1] - r[i]) / (d**(2*j) - 1)
            if i == order - j - 1 and np.abs(r[i] - f_deriv_a(x)) > np.abs(r[i+1] - f_deriv_a(x)):
                return r[i+1]
    return r[0]

y_deriv_rid = np.zeros_like(x)
for i, x_val in enumerate(x):
    y_deriv_rid[i] = ridder(x_val, f, h=.2, d = 2, order = 5)
plt.plot(x, y_deriv_a - y_deriv_rid, label='ridder')
# plt.plot(x, y_deriv_a, label = 'derivative')
# plt.title('functions')
# plt.legend()
# plt.show()

print('ridder:', y_deriv_rid)
print('true:', y_deriv_a)
# plot residuals

# plt.plot(x, y_deriv_a-y_deriv_a, label='true')
plt.title('residuals')
plt.yscale('symlog', linthresh = 1e-16)
plt.legend()
plt.show()

