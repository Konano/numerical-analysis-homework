import numpy as np

eps = 1e-8

def zeroin(f, a, b):
    a, b = np.float64(a), np.float64(b)
    F = lambda x: f(x)
    fa, fb = F(a), F(b)
    if np.sign(fa) == np.sign(fb):
        raise Exception('same sign')
    c, fc = a, fa
    e = d = b - c
    step = 0
    while fb != 0:
        if np.sign(fa) == np.sign(fb):
            a, fa = c, fc
            e = d = b - c
        if np.abs(fa) < np.abs(fb):
            a, b = b, a
            fa, fb = fb, fa
        m = (a - b) / 2
        tol = 2. * eps * max(np.abs(b), 1.)

        if np.abs(m) <= tol or fb == 0.0: # quit
            break
        if np.abs(e) < tol or np.abs(fc) <= np.abs(fb):
            d = e = m
        else:
            s = fb / fc
            if (a == c):
                p = 2. * m * s
                q = 1. - s
            else:
                q, r = fc / fa, fb / fa
                p = s * (2. * m * q * (q - r) - (b - c) * (r - 1.))
                q = (q - 1.) * (r - 1.) * (s - 1.)
            if p > 0:
                q = -q
            else:
                p = -p
            if 2. * p < 3. * m * q - np.abs(tol * q) and p < np.abs(e * q / 2):
                e, d = d, p / q
            else:
                d = e = m
        
        # next iteration
        step += 1
        c, fc = b, fb
        if np.abs(d) > tol:
            b = b + d
        else:
            b = b - np.sign(b - a) * tol
        fb = F(b)
    
    b = np.float128(b)
    print(f'zeroin method took {step} steps to solve the equation: {b:.8f}')
    return b

from mpmath import besselj
import matplotlib.pyplot as plt

j0 = lambda x: besselj(0, x)
x = [.1 * i for i in range(351)]
y = list(map(j0, x))

_, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim((0, 35))
ax.set_ylim((-0.5, 0.5))
ax.set_ylabel('Bessel')
plt.plot(x, y, zorder=2)
plt.grid(True)
plt.axhline(0, color='black', zorder=1)
plt.show()

intervals = [
    (0, 5), (5, 8), (8, 10), (10, 13), (13, 16), 
    (16, 20), (20, 22), (22, 25), (25, 30), (30, 32)
]

zeros = [zeroin(j0, *interval) for interval in intervals]

_, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim((0, 35))
ax.set_ylim((-0.5, 0.5))
ax.set_ylabel('Bessel')
plt.plot(x, y, zorder=2)
plt.grid(True)
plt.axhline(0, color='black', zorder=1)
for zero in zeros:
    plt.scatter(zero, 0, s=20, zorder=3)
plt.savefig('2-3.png')
plt.show()