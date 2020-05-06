import numpy as np
from mpmath import diff, diffs

x = np.array([0.52 , 3.1, 8    , 17.95, 28.65, 39.62, 50.65, 78, 104.6, 156.6, 208.6, 260.7, 312.5, 364.4, 416.3, 468, 494, 507, 520])
y = np.array([5.288, 9.4, 13.84, 20.2 , 24.9 , 28.44, 31.1 , 35, 36.9 , 36.6 , 34.6 , 31   , 26.34, 20.9 , 14.8 , 7.8, 3.7, 1.5, 0.2])

y1_0 =  1.865480
y1_n = -0.046115

sample_x = [2, 30, 130, 350, 515]

h = x[1:] - x[:-1]

def solved(a, b, c, f):
    n = f.shape[0]
    for i in range(1, n):
        b[i] -= a[i-1] / b[i-1] * c[i-1]
        f[i] -= a[i-1] / b[i-1] * f[i-1]
    x = np.zeros_like(f)
    x[-1] = f[-1] / b[-1]
    for i in range(-2, -n-1, -1):
        x[i] = (f[i] - c[i+1] * x[i+1]) / b[i]
    return x

def calM():
    
    lam = np.append([1], h[:-1] / (h[:-1] + h[1:]))
    miu = np.append(h[1:] / (h[:-1] + h[1:]), [1])
    d = np.concatenate((
        [6 / h[0] * ((y[1] - y[0]) / h[0] - y1_0)], 
        6 * (y[:-2] / h[:-1] / (h[:-1] + h[1:]) + y[2:] / h[1:] / (h[:-1] + h[1:]) - y[1:-1] / h[:-1] / h[1:]), 
        [6 / h[-1] * (y1_n - (y[-1] - y[-2]) / h[-1])], 
    ))
    M = solved(miu, np.ones_like(d)*2, lam, d)
    return M

def getS(_x, M):
    i = np.argmax(x > _x) - 1
    return lambda _x: ( \
        M[i] * (x[i+1] - _x)**3 / 6 + \
        M[i+1] * (_x - x[i])**3 / 6 + \
        (y[i] - M[i] * h[i]**2 / 6) * (x[i+1] - _x) + \
        (y[i+1] - M[i+1] * h[i]**2 / 6) * (_x - x[i]) \
    ) / h[i]

M = calM()
S = np.array([list(diffs(getS(x, M), x, 2)) for x in sample_x], dtype=float)
print(S)

import matplotlib.pyplot as plt

___x = np.arange(0.52, 520, 0.005)
___y = np.array([getS(x, M)(x) for x in ___x], dtype=float)
plt.plot(___x, ___y, zorder=1)
plt.scatter(x, y, s=20, zorder=3)
plt.grid(True)
plt.axhline(0, color='black', zorder=1)
plt.savefig('6-8.png')
plt.show()

___x = np.arange(0.52, 520, 0.005)
___y = np.array([diff(getS(x, M), x, 2) for x in ___x], dtype=float)
plt.scatter(___x, ___y, s=1)
plt.grid(True)
plt.axhline(0, color='black', zorder=1)
plt.savefig('6-8-2.png')
plt.show()
