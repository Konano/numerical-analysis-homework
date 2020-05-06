from scipy.linalg import hilbert
import numpy as np

def cholesky(A, size):

    B = np.zeros_like(A)
    for i in range(size):
        B[i, i] = (A[i, i] - ((B[i, :i] ** 2).sum() if i > 0 else 0)) ** .5
        for j in range(i + 1, size):
            B[j, i] = (A[j, i] - B[j, :i].dot(B[i, :i])) / B[i, i]
    return B

def solve(size, delta=0):

    print(f'size = {size}, delta = {delta}')
    H = hilbert(size)
    x = np.zeros(size) + delta + 1
    b = H.dot(x)
    _H = cholesky(H, size)
    xp = np.linalg.solve(_H.T, np.linalg.solve(_H, b))
    print(f'||r||_inf: {abs(b - H.dot(xp)).max()}')
    print(f'||delta||_inf: {abs(xp - x).max()}')

solve(10)
solve(10, 1e-7)
solve(8)
solve(12)