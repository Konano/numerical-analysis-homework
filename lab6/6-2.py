import numpy as np
import matplotlib.pyplot as plt

def solve(R, b, size):
    x = np.zeros_like(b)
    for i in range(size-1, -1, -1):
        for j in range(i+1, size):
            b[i] -= x[j] * R[i, j]
        x[i] = b[i] / R[i, i]
    return x

def Householder(A, f, size):
    for k in range(size):
        sigma = np.sign(A[k, k]) * np.linalg.norm(A[k:, k])
        if sigma == A[k, k]:
            continue
        v = np.zeros(A.shape[0])
        v[k:] = A[k:, k]
        v[k] += sigma
        beta = v.T.dot(v)
        gamma = v.T.dot(f)
        f -= 2 * gamma / beta * v
        for j in range(k, size):
            gamma = v.T.dot(A[:, j])
            A[:, j] -= 2 * gamma / beta * v
    return A, f


def polyval(x, y, D):

    size = len(x)
    A = np.zeros((size, D))
    f = np.zeros((size))
    for i in range(size):
        for j in range(D):
            A[i, j] = x[i] ** j
        f[i] = y[i]

    A, f = Householder(A, f, D)
    x = solve(A[:D], f[:D], D)

    return x

x = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8])
y = np.array([33.4, 79.5, 122.65, 159.05, 214.15, 189.15, 238.65, 252.2, 267.55, 280.5, 296.65, 301.65, 310.4, 318.15, 325.15])

poly_ret = polyval(x, y, 3)
exp_ret = polyval(x, np.log(y), 2)

print('poly:')
print('  a:', poly_ret[0])
print('  b:', poly_ret[1])
print('  c:', poly_ret[2])
print('  R^2:', np.corrcoef(poly_ret[0] + poly_ret[1]*x**1 + poly_ret[2]*x**2, y)[0, 1]**2)
print('exp:')
print('  a:', np.exp(exp_ret[0]))
print('  b:', exp_ret[1])
print('  R^2:', np.corrcoef(np.exp(exp_ret[0] + exp_ret[1] * x), y)[0, 1]**2)