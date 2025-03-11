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

# print('ridder:', y_deriv_rid)
# print('true:', y_deriv_a)
# plot residuals

# plt.plot(x, y_deriv_a-y_deriv_a, label='true')
plt.title('residuals')
plt.yscale('symlog', linthresh = 1e-16)
plt.legend()
# plt.show()
plt.close()

# now you gotta make random numbers on a sphere.
N = 5000
Pu1 = np.random.uniform(0, 1, N)
Pu2 = np.random.uniform(0, 1, N)

theta1 = np.pi * Pu1
phi1 = 2*np.pi * Pu2

theta2 = np.pi*(1 - 2*Pu1)
phi2 = 2*np.pi * Pu2

def to_cartesian(r, theta, phi):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

x1, y1, z1 = to_cartesian(1, theta1, phi1)
x2, y2, z2 = to_cartesian(1, theta2, phi2)
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1, y1, z1, label = 'wrong')
ax.plot(x2, y2, z2, label = 'right')
ax.set_title('random points on a sphere')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.legend()


# random number generator

def xor_shift(seed, a1 = 21, a2 = 35, a3 = 4):
    '''Performs 64bit xor shift on the seed number.'''
    state = np.uint64(seed)
    a1 = np.uint64(a1)
    a2 = np.uint64(a2)
    a3 = np.uint64(a3)

    state = state ^ (state >> a1)
    state = state ^ (state << a2)
    state = state ^ (state >> a3)
    return state

def mult_w_carry(seed, a = 4294957665):
    '''Performs multiply with carry, base 2**32'''
    state = np.uint64(seed)
    a = np.uint64(a)

    state = a*(state & np.uint64(2**32 - 1)) + (state >> np.uint64(32))
    number = state & np.uint64(2**32 - 1)

    return number, state

def rng_int(seed, size):
    state = np.uint64(seed)
    result = np.zeros(size)
    for i in range(size):
        xor_state = xor_shift(state)
        if xor_state < 2**32:
            number, mult_state = mult_w_carry(xor_state)
        else:
            number, mult_state = mult_w_carry(xor_state & np.uint64(2**32 - 1))
        state = mult_state
        result[i] = number

    return result

def rng_float(seed, size, a, b):
    ints = rng_int(seed, size)
    floats = a + (b-a)*(ints / (2**32))
    return floats
plt.close()  

rnds = rng_float(69, 100000, 0, 1)


bins = np.linspace(0, np.max(rnds), 100)
plt.hist(rnds, bins)
plt.show()