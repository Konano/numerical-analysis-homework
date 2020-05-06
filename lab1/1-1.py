import numpy as np
import matplotlib.pyplot as plt

M = 1
eps = 1e-16

def truncating_error(h):
    return M * h / 2

def rounding_error(h):
    return 2 * eps / h

def total_error(h):
    return truncating_error(h) + rounding_error(h)

def actual_error(h):
    return np.abs((np.sin(1 + h) - np.sin(1)) / h - np.cos(1))

x = [10 ** (-16 + i * 0.05) for i in range(321)]

t_errors = list(map(truncating_error, x))
r_errors = list(map(rounding_error, x))
tot_errors = list(map(total_error, x))
act_errors = list(map(actual_error, x))

_, ax = plt.subplots(figsize=(8, 6))
ax.set_xscale("log")
ax.set_xlim((1e-16, 1e0))
ax.set_xlabel('Step h')
ax.set_yscale("log")
ax.set_ylim((1e-17, 1e1))
ax.set_ylabel('Error')
plt.plot(x, t_errors, 'y--', label='Truncating Error')
plt.plot(x, r_errors, 'g--', label='Rounding Error')
plt.plot(x, tot_errors, 'r--', label='Total Error')
plt.plot(x, act_errors, 'b', label='Actual Error')
plt.legend()
plt.savefig('1-1.png')
plt.show()